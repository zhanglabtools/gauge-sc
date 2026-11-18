# Development

The [GAUGE repository](https://github.com/zhanglabtools/gauge-sc) hosts the software, documentation sources and tutorials. See [data release preparation](data-release.md) for versioned tutorial assets.

## Build the documentation locally

From the software repository root:

```bash
python -m pip install -e ".[dev,tutorial,docs]"
python scripts/build_docs.py
python -m http.server 8000 --directory artifacts/docs/html
```

Open `http://localhost:8000`. The build uses Sphinx + Furo with MyST-NB, executes the small-example notebooks in temporary working directories, and treats Sphinx warnings as build failures.

Markdown pages and notebooks remain the single sources of content. The build stages them under ignored `artifacts/docs/source`, preserving relative links. Real-data tutorials retain locally executed outputs with source hashes; use `scripts/check_real_notebooks.py` and the documented data assets to refresh them. A build fails if their notebook sources changed without execution.

Tutorial pages include plots, tables, collapsible setup cells and original notebook downloads. Generated HTML, executed notebooks and datasets do not enter the software distribution.

Read the Docs builds the documentation from this repository using `.readthedocs.yaml`.

## Check the package

```bash
python -m pytest
python scripts/check_notebook.py
python -m build
```

The documentation CI job also executes tutorials. The separate notebook check verifies their standalone execution outside Sphinx.

## Further reading

```{toctree}
:maxdepth: 1

applications
data-release
```
