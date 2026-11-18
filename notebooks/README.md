# Small API examples

For biological applications, start with the [real-data tutorials](../docs/real-data-tutorials.md). The examples below use synthetic data to introduce the concept and individual API operations.

Install with `python -m pip install -e ".[tutorial]"` from the repository root. Each notebook runs independently, creates only small synthetic data in memory, and needs no downloads. Notebook filenames are stable identifiers, not a required reading order.

## Start with the concept

[What does a grouping tell us?](00_concepts.ipynb) follows one example through grouping, information-weighted prediction, complexity cost and view dependence. Read the outputs first; setup and drawing code are marked as collapsible where the viewer supports it.

## Main applications

| Your question | Notebook | Method |
|---|---|---|
| Which clustering granularity should I keep? | [Select a grouping](00_select_grouping.ipynb) | ITCA → ITCA+ |
| What supports different annotations or clusterings? | [Compare fixed groupings](03_compare_groupings.ipynb) | ITCA across views |
| Which groups have relatively greater QC support? | [View-specific audit](01_view_specific_audit.ipynb) | ITCA contribution contrast |
| Are transferred labels supported in the target? | [Evaluate transferred labels](04_transferred_labels.ipynb) | Target-data ITCA (cITCA provenance) |
| Which branches justify finer labels in a spatial dataset? | [Hierarchy and toy space](06_hierarchy_spatial.ipynb) | Conditional ITCA → ancestor-constrained ITCA+ |

## Optional depth

[Screening sensitivity](02_audit_sensitivity.ipynb) reuses one evaluation to compare mirror-ratio settings. [Controlled stress checks](05_stress_checks.ipynb) explores label permutation, cell support and evaluation seeds.

Parameters and return values are in the [API guide](../docs/api.md).

Developer check: `python scripts/check_notebook.py` runs all notebooks in fresh kernels and saves executed copies under ignored `artifacts/`. Tutorials themselves do not export figures or datasets automatically.

[An audit built with your agent](07_agent_audit.ipynb) shows a complete skill-guided workflow: input inspection, analysis choices, actual results and saved-result reuse.
