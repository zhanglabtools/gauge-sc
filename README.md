<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/_static/brand/gauge_sc-logo-dark.svg">
    <source media="(prefers-color-scheme: light)" srcset="docs/_static/brand/gauge_sc-logo-light.svg">
    <img alt="GAUGE" src="docs/_static/brand/gauge_sc-logo-light.svg" width="420">
  </picture>
</p>

GAUGE evaluates biological groupings in single-cell and spatial transcriptomics.
Use ITCA to measure how well a data view predicts group assignments and to
compare molecular and quality-control (QC) contributions. Use ITCA+ to select
a grouping from flat or hierarchical candidates.

- **[Documentation](https://gauge-sc.readthedocs.io/en/latest/)** · Browse the user guide and API reference.
- **[Tutorials](docs/real-data-tutorials.md)** · Learn GAUGE with runnable notebooks and real datasets.
- **[Paper reproduction](reproduction/README.md)** · Replay 41 quantitative panels from frozen results, grouped as Figures 1–6.

## Install

Use Python 3.10 or newer. In an activated virtual environment, install the
alpha release from [PyPI](https://pypi.org/project/gauge-sc/0.1.0a2/) with
either pip or uv.

```bash
python -m pip install gauge-sc==0.1.0a2
```

```bash
uv pip install gauge-sc==0.1.0a2
```

The Python distribution is `gauge-sc`, the module is `gauge_sc`, and the CLI
command is `gauge`. From a local checkout, use `python -m pip install .` or
`python -m pip install ".[tutorial]"` to include the notebook dependencies.
Real-data inputs are available from the [tutorial data release](https://github.com/zhanglabtools/gauge-sc/releases/tag/tutorial-data-v1).

## Quick start

Compare how much the molecular and QC views support each group in a synthetic
dataset that ships with GAUGE.

```python
import gauge_sc as gauge

adata = gauge.datasets.audit_example()
evaluation = gauge.evaluate(
    adata,
    groupby="group",
    views={
        "molecular": {"obsm": "X_pca"},
        "qc": {
            "obs": ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
            "log1p": ["total_counts", "n_genes_by_counts"],
        },
    },
)
audit = gauge.audit(evaluation, reference="molecular", contrast="qc")
print(audit.table)
```

Run the [view-specific audit notebook](notebooks/01_view_specific_audit.ipynb)
for plots, or the [grouping selection notebook](notebooks/00_select_grouping.ipynb)
to try ITCA+.

[Python API](docs/api.md) · [CLI](docs/cli.md) · [Agent skill](docs/agents.md)

## License

GAUGE software is available under the [MIT License](LICENSE). The bundled Inter
fonts use the [SIL Open Font License 1.1](src/gauge_sc/fonts/LICENSE.txt).
