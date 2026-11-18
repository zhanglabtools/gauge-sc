"""Execute all tutorial notebooks in fresh kernels and isolated folders."""

from pathlib import Path
import tempfile

import nbformat
from nbclient import NotebookClient

root = Path(__file__).resolve().parents[1]
for source in sorted((root / "notebooks").glob("*.ipynb")):
    notebook = nbformat.read(source, as_version=4)
    with tempfile.TemporaryDirectory(prefix="gauge_sc-notebook-") as directory:
        NotebookClient(
            notebook,
            timeout=600,
            kernel_name="python3",
            resources={"metadata": {"path": directory}},
        ).execute()
    output = root / "artifacts" / f"{source.stem}.executed.ipynb"
    output.parent.mkdir(exist_ok=True)
    nbformat.write(notebook, output)
    print(f"Executed notebook: {output}")
