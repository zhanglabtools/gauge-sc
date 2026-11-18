"""Version-pinned tutorial assets; no network or disk writes until requested."""

import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
from urllib.request import urlopen


def _catalog():
    return json.loads(Path(__file__).with_name("tutorial_data.json").read_text())


def _cache_root():
    if os.environ.get("GAUGE_CACHE_DIR"):
        return Path(os.environ["GAUGE_CACHE_DIR"]).expanduser()
    if sys.platform == "win32":
        root = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
        return root / "gauge_sc" / "Cache"
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Caches" / "gauge_sc"
    return Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "gauge_sc"


def _verify(path, record):
    if path.stat().st_size != record["bytes"]:
        raise ValueError(
            f"DATA_INTEGRITY: wrong size for {path}. Replace it with the pinned tutorial asset."
        )
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest() != record["sha256"]:
        raise ValueError(
            f"DATA_INTEGRITY: SHA-256 mismatch for {path}. Replace it with the pinned tutorial asset."
        )


def _download(url, target, record):
    """Stage and validate a bounded download before making it visible in cache."""
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=target.parent, suffix=".part", delete=False
        ) as stream:
            temporary = Path(stream.name)
            with urlopen(url, timeout=30) as response:
                received = 0
                while block := response.read(1024 * 1024):
                    received += len(block)
                    if received > record["bytes"]:
                        raise ValueError(
                            "DATA_INTEGRITY: download exceeds expected size."
                        )
                    stream.write(block)
        _verify(temporary, record)
        os.replace(temporary, target)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def load_example(name, *, data_dir=None, cache_dir=None, download=True):
    """Load one checksum-verified, version-pinned real tutorial dataset.

    Args:
        name (str): One of ding_pbmc, ding_mixture, allen_cortex, abc_hpf,
            or vasculature. Only the requested dataset is loaded/downloaded.
        data_dir (path-like, optional): Offline directory holding the release
            h5ad assets. Defaults to GAUGE_TUTORIAL_DATA when set. Explicit
            local input never falls back to a download or trusts a local manifest.
        cache_dir (path-like, optional): Cache root; assets are stored in a
            version subdirectory. Defaults to GAUGE_CACHE_DIR, otherwise
            the platform's user cache location.
        download (bool): Permit download when the pinned cache file is absent.
            False requires an existing cache or local input. Install and import
            never download data. Unpublished assets have no download URL and
            raise a clear error rather than contacting a guessed repository.

    Returns:
        AnnData: An independent in-memory object. The cached file is not changed.
            Expression scale and sampling provenance are retained in uns.

    Raises:
        ValueError: Unknown name, invalid download flag, or corrupt asset.
        FileNotFoundError: Requested local/offline asset is absent.
        RuntimeError: A download was requested before its URL was configured.

    Examples:
        >>> import gauge_sc as gauge
        >>> adata = gauge.datasets.load_example(
        ...     "allen_cortex", data_dir="tutorial-data", download=False
        ... )
    """
    import anndata as ad

    catalog = _catalog()
    if not isinstance(name, str) or name not in catalog["datasets"]:
        raise ValueError(
            f"Unknown tutorial dataset {name!r}; choose from {', '.join(catalog['datasets'])}."
        )
    if not isinstance(download, bool):
        raise ValueError("download must be True or False.")
    record = catalog["datasets"][name]
    local = data_dir if data_dir is not None else os.environ.get("GAUGE_TUTORIAL_DATA")
    if local is not None:
        path = Path(local).expanduser() / record["filename"]
        if not path.is_file():
            raise FileNotFoundError(
                f"Local tutorial asset not found: {path}. No download attempted."
            )
    else:
        root = _cache_root() if cache_dir is None else Path(cache_dir).expanduser()
        path = root / catalog["version"] / record["filename"]
        if not path.is_file():
            if not download:
                raise FileNotFoundError(
                    f"Tutorial asset is not cached: {path}. download=False prevents network access."
                )
            if not record.get("url"):
                raise RuntimeError(
                    "DATA_NOT_PUBLISHED: this tutorial asset has no public URL yet. Supply data_dir or set GAUGE_TUTORIAL_DATA to the prepared release assets."
                )
            _download(record["url"], path, record)
    _verify(path, record)
    data = ad.read_h5ad(path)
    if data.shape != (record["n_obs"], record["n_vars"]):
        raise ValueError("DATA_INTEGRITY: asset dimensions disagree with the catalog.")
    return data
