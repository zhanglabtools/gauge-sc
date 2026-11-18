# Applications, APIs and tutorials

Start with the [concept guide](concepts.md), then choose a [real-data tutorial](real-data-tutorials.md) for your question. The small examples below focus on individual API operations.

ITCA measures prediction support; ITCA+ balances that support against grouping complexity. Reference, inferred and transferred labels use the same evaluation API.

| Application | API | Tutorial / status | Validation |
|---|---|---|---|
| Select clustering granularity | `evaluate` → `select` | [Small synthetic selection](../notebooks/00_select_grouping.ipynb), implemented for fixed flat candidates | Formula, decision transitions, monotone complexity, identical cells, deterministic ties |
| View-specific QC audit | `evaluate` → `audit` | [Introductory audit](../notebooks/01_view_specific_audit.ipynb), implemented | Paired views, contributions, screening reasons and result export |
| Screening parameter sensitivity | Reuse one evaluation with `audit` | [Advanced example](../notebooks/02_audit_sensitivity.ipynb), implemented | No refitting, fixed scores, candidate sets |
| Compare reference and inferred labels across views | `evaluate` | [Compare fixed groupings](../notebooks/03_compare_groupings.ipynb), implemented | Fixed candidate provenance, explicit views and common observations |
| Evaluate transferred labels in target data | `evaluate` on target AnnData | [Target-label evaluation](../notebooks/04_transferred_labels.ipynb), implemented; transfer remains external | Target proportions and target CV; no source/target score conflation |
| Select hierarchical / spatial structures | `evaluate(hierarchy=...)` → `select` → explicit writeback | [Hierarchy and toy space](../notebooks/06_hierarchy_spatial.ipynb), implemented | Exhaustive small-tree optimum, research-DP parity, unchanged global assignments; toy space is not anatomical validation |
| Compare workflow candidates | `evaluate` → `select` on common observations and a shared view | [Fixed workflow endpoints](../notebooks/real_data/13_workflow_endpoints.ipynb) | Apply the same evaluation protocol to each candidate and report small-group support separately from selection |
| Stress tests | Evaluation under controlled perturbations | [Controlled stress checks](../notebooks/05_stress_checks.ipynb), implemented; broader classifier tests remain optional | Fixed proportions under label permutation, explicit CV failure, nested subsets and seed sensitivity; no full benchmark grid |

For parameters and return values, use the [API guide](api.md) and [Python reference](reference.rst).
