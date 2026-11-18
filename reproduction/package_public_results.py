"""Build a portable quantitative reproduction asset from a verified local freeze."""

import argparse
from collections import defaultdict
import csv
from importlib.metadata import version
import json
from pathlib import Path
import shutil
import tempfile

import numpy as np
import pandas as pd

from render import verify
from source_asset import (
    COSMX_ACCESS,
    COSMX_ARCHIVE,
    COSMX_METADATA,
    COSMX_PANELS,
    PANEL_COUNTS,
    digest,
)

SCRIPTS = ["render.py", "validate.py", "source_asset.py", "fetch_cosmx_metadata.py"]
DEPENDENCIES = [
    "numpy",
    "pandas",
    "matplotlib",
    "seaborn",
    "scipy",
    "scikit-learn",
    "requests",
    "pillow",
]
PATH_COLUMNS = {"source_table", "source_table_path", "source_h5ad", "hvg_list_source"}
ID_PANELS = {"3", "4", "5", "6", "7"}
PLOT_COLUMNS = [
    "left_plot_category",
    "left_plot_color",
    "right_plot_category",
    "right_plot_color",
]


def copy_result_table(source, output, identifiers):
    with (
        source.open(newline="", encoding="utf-8") as original,
        output.open("w", newline="", encoding="utf-8") as target,
    ):
        reader = csv.DictReader(original)
        columns = [name for name in reader.fieldnames if name not in PATH_COLUMNS]
        writer = csv.DictWriter(target, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            if identifiers is not None and "cell_id" in columns:
                row["cell_id"] = identifiers.setdefault(
                    row["cell_id"], f"cell{len(identifiers) + 1}"
                )
            writer.writerow(row)


def copy_cosmx_results(source, output, metadata):
    table = pd.read_csv(source)
    selected = metadata.loc[table.cell_id]
    np.testing.assert_array_equal(
        selected[["CenterX_global_px", "CenterY_global_px"]], table[["x", "y"]]
    )
    results = table[PLOT_COLUMNS].copy()
    results.insert(0, "metadata_row", selected.metadata_row.to_numpy())
    results.to_csv(output, index=False)


def package(bundle, cosmx_metadata, output):
    manifest = verify(bundle)
    if (
        manifest["schema_version"] != 2
        or manifest.get("figures", list(PANEL_COUNTS)) != list(PANEL_COUNTS)
        or manifest.get("external_inputs")
    ):
        raise ValueError("Use the full local quantitative asset.")
    if output.exists():
        raise ValueError("Use a new output directory.")

    metadata = pd.read_csv(
        cosmx_metadata,
        usecols=["fov", "cell_ID", "CenterX_global_px", "CenterY_global_px"],
    )
    metadata["metadata_row"] = metadata.index
    metadata.index = metadata.fov.astype(str) + "_" + metadata.cell_ID.astype(str)
    if not metadata.index.is_unique:
        raise ValueError("CosMx metadata has duplicate cell identifiers.")

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as directory:
        root = Path(directory) / "asset"
        source_data = root / "source_data"
        source_data.mkdir(parents=True)
        reference = root / "reference"
        reference.mkdir()
        scripts = root / "reproduction"
        scripts.mkdir()
        identifiers = defaultdict(dict)
        names = sorted(
            {name for panel in manifest["panels"] for name in panel["tables"]}
        )
        for name in names:
            source = bundle / "source_data" / name
            dest = source_data / name
            if source.stem in COSMX_PANELS:
                copy_cosmx_results(source, dest, metadata)
            else:
                ids = identifiers[name[3]] if name[3] in ID_PANELS else None
                copy_result_table(source, dest, ids)
        shutil.copyfile(
            bundle / "reference/fig2h_per_dataset.csv",
            reference / "fig2h_per_dataset.csv",
        )
        here = Path(__file__).resolve().parent
        for name in SCRIPTS:
            shutil.copyfile(here / name, scripts / name)
        shutil.copyfile(here / "PUBLIC-README.md", root / "README.md")
        shutil.copyfile(here / "DATA-SOURCES.md", root / "DATA-SOURCES.md")
        shutil.copyfile(here.parent / "LICENSE", root / "LICENSE")
        (root / "requirements.txt").write_text(
            "".join(f"{name}=={version(name)}\n" for name in DEPENDENCIES)
        )
        public = {
            "schema_version": 2,
            "scope": "Frozen source-data reconstruction of 41 quantitative panels",
            "panel_numbering": "paper",
            "panels": manifest["panels"],
            "external_inputs": {
                "cosmx_metadata": {
                    "panels": sorted(COSMX_PANELS),
                    "filename": COSMX_METADATA,
                    "bytes": cosmx_metadata.stat().st_size,
                    "sha256": digest(cosmx_metadata),
                    "access_url": COSMX_ACCESS,
                    "archive_url": COSMX_ARCHIVE,
                }
            },
            "files": {
                path.relative_to(root).as_posix(): {
                    "bytes": path.stat().st_size,
                    "sha256": digest(path),
                }
                for path in sorted(root.rglob("*"))
                if path.is_file()
            },
        }
        (root / "manifest.json").write_text(json.dumps(public, indent=2) + "\n")
        verify(root)
        root.rename(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--cosmx-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    package(args.bundle, args.cosmx_metadata, args.output)
    print(f"Packaged 41 quantitative panels at {args.output}")
