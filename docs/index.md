```{image} docs/_static/brand/gauge_sc-logo-light.svg
:class: project-logo only-light
:alt: GAUGE
:width: 360px
```

```{image} docs/_static/brand/gauge_sc-logo-dark.svg
:class: project-logo only-dark
:alt: GAUGE
:width: 360px
```

# GAUGE: which distinctions are worth keeping?

A grouping proposes distinctions between cells. **ITCA measures how much information a specified data view can predict about those distinctions. ITCA+ balances that support against the cost of retaining finer groups.**

Start with one small example: the same cells, three groupings, and two data views.

[**Understand the concept →**](docs/concepts.md) · [Run the concept notebook](notebooks/00_concepts.ipynb)

![The same observations under three grouping structures](docs/assets/concept-groupings.png)

## Apply the concept to real data

```{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} Audit groups for QC effects
<img src="docs/assets/application-audit.svg" class="application-icon" alt="Two-view profiling and group audit">

[Human–mouse mixture](notebooks/real_data/11_mixture_audit.ipynb)
:::
:::{grid-item-card} Select grouping granularity
<img src="docs/assets/application-selection.svg" class="application-icon" alt="Granularity changes with the complexity penalty">

[Allen Cortex selection](notebooks/real_data/14_allen_granularity.ipynb)
:::
:::{grid-item-card} Adapt a spatial taxonomy
<img src="docs/assets/application-hierarchy.svg" class="application-icon" alt="Adaptive taxonomy mapped to spatial observations">

[ABC MERFISH hierarchy](notebooks/real_data/12_spatial_hierarchy.ipynb)
:::
:::{grid-item-card} Advanced comparisons
<img src="docs/assets/application-workflow.svg" class="application-icon" alt="Selection among candidate workflow endpoints">

[PBMC platform profiles](notebooks/real_data/10_platform_profiles.ipynb) · [Brain-vasculature workflow candidates](notebooks/real_data/13_workflow_endpoints.ipynb)
:::
```

See the [real-data tutorial directory and data setup](docs/real-data-tutorials.md). Application icons are adapted from Figure 1c.

## Reproduce the paper figures

The [paper reproduction guide](docs/paper-reproduction.md) explains how to redraw 41 quantitative panels from frozen results. The [standalone Figure 4 replay](docs/reproduce-figure4.md) redraws panels 4b–f and verifies both mixture results from pinned source tables.

## Three computational entry points

**`evaluate` → ITCA support. `audit` → view contrasts. `select` → ITCA+ structure selection.**

Use AnnData inputs, inspect result tables, and display selected labels or audit scores with Scanpy. Parameters and validity rules are in the [API guide](docs/api.md); exact signatures are in the [Python reference](docs/reference.rst).

The concept guide uses a small synthetic example; application tutorials take the same ideas to real data.

```{note}
Get the software and tutorial datasets from the [public GAUGE repository](https://github.com/zhanglabtools/gauge-sc). See [installation](docs/installation.md) and [tutorial data setup](docs/real-data-tutorials.md) to run the notebooks locally.
```

```{toctree}
:hidden:
:caption: Start here

Concept guide <docs/concepts>
Installation <docs/installation>
Work with your agent <docs/agents>
```

```{toctree}
:hidden:
:caption: Real-data applications

Tutorial directory and data setup <docs/real-data-tutorials>
Mixture audit <notebooks/real_data/11_mixture_audit>
Allen Cortex selection <notebooks/real_data/14_allen_granularity>
Spatial hierarchy <notebooks/real_data/12_spatial_hierarchy>
```

```{toctree}
:hidden:
:caption: Advanced real-data comparisons

Platform profiles <notebooks/real_data/10_platform_profiles>
Workflow endpoints <notebooks/real_data/13_workflow_endpoints>
```

```{toctree}
:hidden:
:caption: Small API examples

Choose an example <notebooks/README>
Concept in code <notebooks/00_concepts>
Select a grouping <notebooks/00_select_grouping>
Compare groupings <notebooks/03_compare_groupings>
View-specific audit <notebooks/01_view_specific_audit>
Transferred labels <notebooks/04_transferred_labels>
Hierarchy and toy space <notebooks/06_hierarchy_spatial>
An agent-generated audit <notebooks/07_agent_audit>
```

```{toctree}
:hidden:
:caption: Go deeper

Screening sensitivity <notebooks/02_audit_sensitivity>
Controlled stress checks <notebooks/05_stress_checks>
Recipes <docs/recipes>
API guide <docs/api>
CLI guide <docs/cli>
Python reference <docs/reference>
Data and provenance <docs/data>
Paper reproduction <docs/paper-reproduction>
Numerical reproducibility <docs/numerical-reproducibility>
```

```{toctree}
:hidden:
:caption: Project

Development <docs/development>
```
