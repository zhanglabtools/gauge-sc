"""Extract a self-contained Figure 4 source asset from a verified main-figure asset."""

import argparse
import json
from pathlib import Path
import re
import shutil
import tempfile

from render import verify
from source_asset import SCHEMATICS, digest


def package(bundle, output):
    manifest = verify(bundle)
    if manifest.get("figures") == [4]:
        raise ValueError("Use the full main-figure asset as input.")
    if output.exists():
        raise ValueError("Use a new output directory; frozen evidence is immutable.")

    panels = [
        panel
        for panel in manifest["panels"]
        if panel["figure"] == 4 and panel["id"] not in SCHEMATICS
    ]
    inputs = {f"source_data/{name}" for panel in panels for name in panel["tables"]}
    limitations = []
    for limitation in manifest.get("limitations", []):
        match = re.match(r"Figure (\d+)", limitation)
        if match is None or match[1] == "4":
            limitations.append(limitation)

    scoped_manifest = {
        "schema_version": 2,
        "figures": [4],
        "scope": "Figure 4b–f source-data reconstruction",
        "source_asset_manifest_sha256": digest(bundle / "manifest.json"),
        "panels": panels,
        "files": {name: manifest["files"][name] for name in sorted(inputs)},
    }
    if manifest.get("panel_numbering") == "paper":
        scoped_manifest["panel_numbering"] = "paper"
    else:
        scoped_manifest.update(
            authority=manifest["authority"],
            research_revision=manifest.get("research_revision"),
            research_working_tree_clean=manifest.get("research_working_tree_clean"),
            limitations=limitations,
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as staging:
        staged = Path(staging) / "source"
        for name in scoped_manifest["files"]:
            destination = staged / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(bundle / name, destination)
        (staged / "manifest.json").write_text(
            json.dumps(scoped_manifest, indent=2) + "\n"
        )
        verify(staged)
        staged.rename(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    package(args.bundle, args.output)
    print(f"Packaged Figure 4 source data at {args.output}")
