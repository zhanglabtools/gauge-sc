"""Pinned dataset loading must stay offline unless a cache miss permits download."""

import hashlib
import io

import anndata as ad
import numpy as np
import pytest

from gauge_sc import datasets
from gauge_sc import _tutorial_data as assets


@pytest.fixture
def example(tmp_path, monkeypatch):
    folder = tmp_path / "local"
    folder.mkdir()
    path = folder / "tiny.h5ad"
    ad.AnnData(np.eye(3)).write_h5ad(path)
    payload = path.read_bytes()
    catalog = {
        "version": "tutorial-data-test",
        "datasets": {
            "tiny": {
                "filename": "tiny.h5ad",
                "bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
                "n_obs": 3,
                "n_vars": 3,
                "url": "https://example.invalid/tiny.h5ad",
            }
        },
    }
    monkeypatch.setattr(assets, "_catalog", lambda: catalog)
    monkeypatch.delenv("GAUGE_TUTORIAL_DATA", raising=False)
    return folder, payload, catalog


def no_network(*args, **kwargs):
    raise AssertionError("Unexpected network access")


def test_local_load_is_offline_and_returns_independent_objects(example, monkeypatch):
    folder, _, _ = example
    monkeypatch.setattr(assets, "urlopen", no_network)
    first = datasets.load_example("tiny", data_dir=folder)
    first.X[:] = 8
    np.testing.assert_array_equal(
        datasets.load_example("tiny", data_dir=folder).X, np.eye(3)
    )
    monkeypatch.setenv("GAUGE_TUTORIAL_DATA", str(folder))
    assert datasets.load_example("tiny", download=False).shape == (3, 3)


def test_download_then_versioned_cache_is_reused(example, tmp_path, monkeypatch):
    _, payload, _ = example
    calls = []

    def fetch(url, timeout):
        calls.append(url)
        return io.BytesIO(payload)

    monkeypatch.setattr(assets, "urlopen", fetch)
    cache = tmp_path / "cache"
    assert datasets.load_example("tiny", cache_dir=cache).shape == (3, 3)
    monkeypatch.setattr(assets, "urlopen", no_network)
    assert datasets.load_example("tiny", cache_dir=cache, download=False).shape == (
        3,
        3,
    )
    assert (cache / "tutorial-data-test" / "tiny.h5ad").is_file()
    assert len(calls) == 1


@pytest.mark.parametrize("payload_kind", ["truncated", "oversized", "wrong_hash"])
def test_bad_downloads_do_not_enter_cache(example, tmp_path, monkeypatch, payload_kind):
    _, payload, _ = example
    broken = {
        "truncated": payload[:20],
        "oversized": payload + b"extra",
        "wrong_hash": b"x" * len(payload),
    }[payload_kind]
    monkeypatch.setattr(assets, "urlopen", lambda *args, **kwargs: io.BytesIO(broken))
    cache = tmp_path / "cache"
    with pytest.raises(ValueError, match="DATA_INTEGRITY"):
        datasets.load_example("tiny", cache_dir=cache)
    assert not list(cache.rglob("*.h5ad"))
    assert not list(cache.rglob("*.part"))


def test_missing_local_input_never_falls_back_to_download(
    example, tmp_path, monkeypatch
):
    monkeypatch.setattr(assets, "urlopen", no_network)
    with pytest.raises(FileNotFoundError, match="No download attempted"):
        datasets.load_example("tiny", data_dir=tmp_path / "missing")


def test_corrupt_local_input_is_rejected(example, monkeypatch):
    folder, payload, _ = example
    (folder / "tiny.h5ad").write_bytes(b"x" * len(payload))
    monkeypatch.setattr(assets, "urlopen", no_network)
    with pytest.raises(ValueError, match="SHA-256"):
        datasets.load_example("tiny", data_dir=folder)


def test_offline_cache_miss_and_unpublished_asset_are_explicit(
    example, tmp_path, monkeypatch
):
    _, _, catalog = example
    monkeypatch.setattr(assets, "urlopen", no_network)
    with pytest.raises(FileNotFoundError, match="download=False"):
        datasets.load_example("tiny", cache_dir=tmp_path, download=False)
    catalog["datasets"]["tiny"]["url"] = None
    with pytest.raises(RuntimeError, match="DATA_NOT_PUBLISHED"):
        datasets.load_example("tiny", cache_dir=tmp_path)
    assert not (tmp_path / "tutorial-data-test").exists()


def test_dataset_name_and_download_flag_validation(example, tmp_path):
    with pytest.raises(ValueError, match="Unknown tutorial"):
        datasets.load_example("../tiny", cache_dir=tmp_path)
    with pytest.raises(ValueError, match="download must"):
        datasets.load_example("tiny", cache_dir=tmp_path, download="false")
