"""Explicit developer step: extract concept figures from a checked notebook.

Run scripts/check_notebook.py first. This never executes analysis or exports data.
The tutorial itself does not save files; these small PNGs illustrate the docs page.
"""

import base64
from pathlib import Path

import nbformat

root = Path(__file__).resolve().parents[1]
source = nbformat.read(root / "notebooks/00_concepts.ipynb", as_version=4)
executed = nbformat.read(root / "artifacts/00_concepts.executed.ipynb", as_version=4)
if [cell.source for cell in source.cells] != [cell.source for cell in executed.cells]:
    raise ValueError("Executed concept notebook is stale; run check_notebook.py first.")
output = root / "docs/assets"
output.mkdir(exist_ok=True)
for cell in executed.cells:
    name = cell.metadata.get("concept_figure")
    if not name:
        continue
    images = [
        item["data"]["image/png"]
        for item in cell.get("outputs", [])
        if "image/png" in item.get("data", {})
    ]
    if len(images) != 1:
        raise ValueError(f"Expected exactly one figure for {name!r}.")
    (output / name).write_bytes(base64.b64decode(images[0]))
    print(f"Updated {output / name}")
