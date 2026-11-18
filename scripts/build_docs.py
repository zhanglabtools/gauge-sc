"""Stage single-source tutorials and build the local Sphinx/Furo website."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

import nbformat


def render_details(text):
    """Let MyST parse Markdown inside portable HTML details blocks."""
    return re.sub(
        r"<details>\s*<summary>(.*?)</summary>\s*(.*?)\s*</details>",
        lambda match: ("\n:::{dropdown} " + match[1] + "\n\n" + match[2] + "\n:::\n"),
        text,
        flags=re.DOTALL,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    destination = root / "artifacts" / "docs"
    source = destination / "source"
    if source.exists():
        shutil.rmtree(source)
    source.mkdir(parents=True)
    shutil.copytree(root / "docs", source / "docs")
    shutil.move(source / "docs" / "conf.py", source / "conf.py")
    shutil.move(source / "docs" / "index.md", source / "index.md")
    shutil.copytree(
        root / "notebooks",
        source / "notebooks",
        ignore=shutil.ignore_patterns(".ipynb_checkpoints", "outputs", "__pycache__"),
    )
    downloads = source / "downloads"
    downloads.mkdir()
    for path in source.rglob("*.md"):
        path.write_text(render_details(path.read_text()))
    for path in sorted((source / "notebooks").rglob("*.ipynb")):
        shutil.copy2(root / path.relative_to(source), downloads / path.name)
        notebook = nbformat.read(path, as_version=4)
        if path.parent.name == "real_data":
            actual = hashlib.sha256(
                json.dumps([c.source for c in notebook.cells]).encode()
            ).hexdigest()
            if (
                notebook.metadata.get("gauge_sc_execution", {}).get("source_sha256")
                != actual
            ):
                raise RuntimeError(
                    f"Real-data outputs are missing or stale: {path.name}. Run scripts/check_real_notebooks.py."
                )

        for cell in notebook.cells:
            if cell.cell_type == "markdown":
                cell.source = render_details(cell.source)
        notebook.cells.insert(
            1,
            nbformat.v4.new_markdown_cell(
                "{download}`Download this notebook </downloads/" + path.name + ">`"
            ),
        )
        nbformat.write(notebook, path)
    print(f"Staged documentation: {source}", flush=True)
    if not args.prepare_only:
        subprocess.run(
            [
                sys.executable,
                "-m",
                "sphinx",
                "-b",
                "html",
                "-W",
                "--keep-going",
                str(source),
                str(destination / "html"),
            ],
            cwd=root,
            check=True,
        )
        reference = (destination / "html" / "docs" / "reference.html").read_text()
        for name in (
            "gauge_sc.evaluate",
            "gauge_sc.audit",
            "gauge_sc.select",
            "gauge_sc.pl.audit",
            "gauge_sc.pl.table",
            "gauge_sc.pl.hierarchy",
            "gauge_sc.datasets.audit_example",
            "gauge_sc.datasets.load_example",
            "gauge_sc.io.load_result",
        ):
            if f'id="{name}"' not in reference:
                raise RuntimeError(f"API reference did not render a signature: {name}")
        if ".. py:function::" in reference or ".. py:class::" in reference:
            raise RuntimeError("API reference contains unparsed autodoc directives.")
        print(f"Local website: {destination / 'html' / 'index.html'}")


if __name__ == "__main__":
    main()
