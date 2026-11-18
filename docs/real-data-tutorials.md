# Real-data tutorials

Choose a real-data application below. Each tutorial introduces a question, runs the public API and explains the resulting plots. The [concept notebook](../notebooks/00_concepts.ipynb) provides a short introduction that runs without a dataset download.

Start with the [mixture audit](../notebooks/real_data/11_mixture_audit.ipynb). Its notebook contains the install commands, instructions for loading `ding_mixture.h5ad` beside the notebook, a guide to reading the output, and a runnable example you can adapt to your own AnnData. Follow the [paper reproduction guide](paper-reproduction.md) for panel-specific methods, source data and validation.

Choose [Allen granularity](../notebooks/real_data/14_allen_granularity.ipynb) to compare clustering resolutions on the same cells. Its self-contained notebook includes a two-row score comparison at `δ=0.02` and explains when numerical differences can change the selected granularity.

Choose [spatial hierarchy](../notebooks/real_data/12_spatial_hierarchy.ipynb) to vary the annotation depth across a supplied taxonomy. The notebook uses `abc_hpf.h5ad`, walks through the DG, VLMC and CA3 branch decisions, and maps the selected depths back to the same tissue coordinates.

For advanced comparisons, use [PBMC platform profiles](../notebooks/real_data/10_platform_profiles.ipynb) to read expression and measured-QC support within each platform, or [workflow endpoints](../notebooks/real_data/13_workflow_endpoints.ipynb) to compare three fixed candidate partitions using a common evaluation view. Each notebook includes local data setup and a runnable starting point for your own AnnData.

| Manuscript | Core readout retained | Tutorial | Scope |
|---|---|---|---|
| Figure 2 | A small set of perturbation response curves | Real-data extension pending | Advanced perturbation analysis |
| Figure 3d–e | Molecular and QC ITCA profiles over grouping size | [Platform profiles](../notebooks/real_data/10_platform_profiles.ipynb) | Two PBMC1 platforms, 1,000 cells each; three resolutions |
| [Figure 4b–d](reproduce-figure4.md) | Contrast map, prespecified follow-up, and withheld-species readout | [Mixture audit](../notebooks/real_data/11_mixture_audit.ipynb) | Complete 3,159-cell Mixture 1 10x sample; one Leiden grouping |
| Figure 5d,f and a selected-label map | Raw support, selection path and selected labels | [Allen granularity](../notebooks/real_data/14_allen_granularity.ipynb) | 3,000 uniformly sampled cells; 2,000 HVGs, five Leiden candidates |
| Figure 6b and spatial readouts | Pruning path and real spatial label map | [Spatial hierarchy](../notebooks/real_data/12_spatial_hierarchy.ipynb) | 3,000 HPF slice cells; supplied taxonomy |
| Figure 7a,e | Fixed candidate evaluation and population-support summary | [Workflow endpoints](../notebooks/real_data/13_workflow_endpoints.ipynb) | 2,500 cells; three Scanpy candidates |

## Tutorial datasets

The [tutorial-data-v1 release](https://github.com/zhanglabtools/gauge-sc/releases/tag/tutorial-data-v1) provides five `.h5ad` files. See [data sources and licenses](https://github.com/zhanglabtools/gauge-sc/blob/tutorial-data-v1/DATA-SOURCES.md) for attribution and terms.

| Dataset key | Cells | Features | Download MB |
|---|---:|---:|---:|
| `ding_pbmc` | 2,000 | 33,694 | 6.53 |
| `ding_mixture` | 3,159 | 15,796 | 28.15 |
| `abc_hpf` | 3,000 | 500 | 1.50 |
| `vasculature` | 2,500 | 26,666 | 14.48 |
| `allen_cortex` | 3,000 | 2,000 | 2.19 |

All five assets total about 52.8 MB. They are separate Release attachments, excluded from Git history and Python distributions.

## Load a dataset in Python

`load_example` downloads the requested dataset on first use and reuses the cached copy on later calls. Once the file is cached, pass `download=False` to run offline.

```python
import gauge_sc as gauge

adata = gauge.datasets.load_example("allen_cortex")
adata = gauge.datasets.load_example("allen_cortex", download=False)
```

Use `cache_dir` or `GAUGE_CACHE_DIR` to choose the cache location. Files are stored under its `tutorial-data-v1` subdirectory.

## Run a notebook with a local file

Download the `.h5ad` listed in the notebook and start Jupyter from the folder containing it. For a shared folder, set `GAUGE_TUTORIAL_DATA` to its absolute path or edit `data_dir` in the notebook. You can use the same file in a script by passing its directory explicitly.

```python
import gauge_sc as gauge

adata = gauge.datasets.load_example("allen_cortex", data_dir=".", download=False)
```

## Prepare assets from manuscript inputs

Developers can reproduce the fixed subset recipe without downloading or modifying source data:

```bash
python -m pip install -e ".[tutorial]"
python scripts/prepare_tutorial_data.py --research-root /path/to/research-repository --allen-input /path/to/allen_cortex_20k.h5ad --output artifacts/tutorial-data
python scripts/package_tutorial_data.py --input artifacts/tutorial-data --output artifacts/data-releases/tutorial-data-v1 --version tutorial-data-v1 --github-repo zhanglabtools/gauge-sc
export GAUGE_TUTORIAL_DATA="$PWD/artifacts/tutorial-data"
python scripts/check_real_notebooks.py
```

To refresh one tutorial, pass its notebook stem to `--notebook`, for example `--notebook 12_spatial_hierarchy`. The notebook records the runtime and input manifest used for the saved outputs.

Output directories must be new. The package catalog pins the reviewed release bytes; if regenerated assets differ, review and version them before updating that catalog. The packaging script only copies validated files and writes metadata; it never uploads. See [maintainer release preparation](data-release.md).

Sampling uses seed 0 independently per dataset, without balancing biological labels. Species annotations enter only the post-audit readout. Allen normalization uses the full gene set before HVG restriction, and its tutorial does not normalize or log again. ABC spatial expression is supplied log2 data. Original QC covariates retain their full-gene-set provenance.

Each asset records source identity and checksum, sampling rules, matrix dimensions, byte size and preprocessing. Preserve upstream attribution and redistribution terms when publishing; no new dataset license is assigned by this software.
