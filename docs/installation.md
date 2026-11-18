# Installation

Get GAUGE from the public [zhanglabtools/gauge-sc](https://github.com/zhanglabtools/gauge-sc) repository. Use [GitHub Releases](https://github.com/zhanglabtools/gauge-sc/releases) to find a tagged version or its release assets.

Use Python 3.10 or newer. To install from source, clone the repository and create an isolated environment.

```bash
git clone https://github.com/zhanglabtools/gauge-sc.git
cd gauge-sc
python -m venv .venv
source .venv/bin/activate
python -m pip install ".[tutorial]"
```

The `tutorial` extra includes Scanpy and the notebook execution dependencies. From the source directory, use `python -m pip install .` to install the core API. To install a downloaded wheel, pass its local path to `python -m pip install`.

Verify the installation:

```python
import gauge_sc as gauge

print(gauge.__version__)
adata = gauge.datasets.audit_example()
print(adata.shape)
```

Then start with [the concept guide](concepts.md) and [its runnable notebook](../notebooks/00_concepts.ipynb). Each tutorial page provides a notebook download. Open it in your Jupyter environment with this Python environment selected as its kernel. To install JupyterLab locally, run `python -m pip install jupyterlab`.

## Data and optional dependencies

Installation does not download datasets. Small API examples generate synthetic data; real-data tutorials use separate dataset assets. See [tutorial data setup](real-data-tutorials.md).

The `validation` extra adds optional dependencies for research validation. It is not required for the main tutorials. There is no dependency on the legacy `itca` or `scit` packages.
