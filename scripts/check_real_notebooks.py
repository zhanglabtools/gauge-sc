"""Execute real-data tutorials using explicitly supplied local dataset assets."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
from textwrap import dedent

import nbformat
from nbclient import NotebookClient

RUNTIME_CELL = dedent("""
    from importlib.metadata import version
    import json
    from pathlib import Path
    import platform

    import numpy
    import scipy
    from threadpoolctl import threadpool_info

    threadpools = threadpool_info()
    for library in threadpools:
        library["filepath"] = Path(library["filepath"]).name
    build_blas = {}
    for package in (numpy, scipy):
        config = getattr(package.__config__, "CONFIG", {})
        blas = config.get("Build Dependencies", {}).get("blas", {})
        build_blas[package.__name__] = {
            key: blas[key] for key in ("name", "version") if key in blas
        }
    packages = ["gauge_sc", "scanpy", "anndata", "scikit-learn", "numpy", "scipy"]
    print(json.dumps({
        "python": platform.python_version(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "threadpools": threadpools,
        "build_blas": build_blas,
        "versions": {name: version(name) for name in packages},
    }))
""").strip()


def main():
    root = Path(__file__).resolve().parents[1]
    sources = {
        p.stem: p for p in sorted((root / "notebooks" / "real_data").glob("*.ipynb"))
    }
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--notebook",
        choices=sources,
        action="append",
        help="Notebook stem to execute. Defaults to all real-data notebooks.",
    )
    args = parser.parse_args()
    if "GAUGE_TUTORIAL_DATA" not in os.environ:
        raise SystemExit(
            "Set GAUGE_TUTORIAL_DATA to the prepared asset directory; see docs/real-data-tutorials.md."
        )
    data = Path(os.environ["GAUGE_TUTORIAL_DATA"]).resolve()
    os.environ["GAUGE_TUTORIAL_DATA"] = str(data)
    manifest = data / "manifest.json"
    for name in dict.fromkeys(args.notebook or sources):
        source = sources[name]
        notebook = nbformat.read(source, as_version=4)
        notebook.cells.append(nbformat.v4.new_code_cell(RUNTIME_CELL))
        with tempfile.TemporaryDirectory(prefix="gauge_sc-real-notebook-") as directory:
            NotebookClient(
                notebook,
                timeout=600,
                kernel_name="python3",
                resources={"metadata": {"path": directory}},
            ).execute()
        runtime_cell = notebook.cells.pop()
        runtime = json.loads(
            "".join(
                output.get("text", "")
                for output in runtime_cell.outputs
                if output.output_type == "stream" and output.name == "stdout"
            )
        )
        notebook.metadata["mystnb"] = {"execution_mode": "off"}
        notebook.metadata["gauge_sc_execution"] = {
            **runtime,
            "data_manifest_sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
            "source_sha256": hashlib.sha256(
                json.dumps([c.source for c in notebook.cells]).encode()
            ).hexdigest(),
        }
        # Only compact rendered outputs enter Git; input h5ad files never do.
        nbformat.write(notebook, source)
        print(f"Executed real-data tutorial: {source.name}", flush=True)


if __name__ == "__main__":
    main()
