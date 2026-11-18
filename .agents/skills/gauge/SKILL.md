---
name: gauge
description: "Use ITCA+ to evaluate prediction support for supplied biological groupings, audit group contributions across data views, or select flat or hierarchical grouping granularity in single-cell and spatial AnnData. Use for ITCA/ITCA+ analysis requests, interpretation of saved results, and Scanpy integration."
---

# GAUGE

ITCA measures information-weighted predictability of supplied group distinctions
from a specified view. ITCA+ balances that support against retained complexity.
Choose an analysis that answers the user's question; this skill does not fix a
classifier, view, preprocessing recipe, screening ratio, or complexity penalty.

## Frame the task

- **Evaluate**: compare support for fixed groupings, labels or workflow endpoints.
- **Audit**: contrast group contributions between two views for one grouping.
- **Select**: choose among supplied candidates, or prune a supplied hierarchy.

Infer intent from the request and available tutorial/context. Inspect the data
before choosing settings. If the meaning of labels, the reference view, or the
selection goal is unclear, ask a concise scientific question. Make routine
execution choices yourself and explain consequential analysis choices.

## Inspect and decide

1. Check the active Python environment, installed ITCA+ version and public API.
   Skill installation does not install the Python package. Use the existing
   analysis environment or install the user's software checkout when needed.
2. Inventory `.obs`, `.obsm`, layers, cell identities, group support and preprocessing
   provenance. The CLI `inspect` command helps; inspect values in Python as needed.
   Do not infer that `.X` is normalized, that PCA is suitable, or that QC columns
   exist from their names alone. Reuse a saved evaluation if it matches the task.
3. Choose the grouping(s) and named views with a stated rationale. Generating
   clusters, embeddings or workflow candidates is a separate Scanpy step.
   Preserve a common observation set when comparing candidates. Review control
   views separately and choose a selection view explicitly.
4. Select evaluation settings appropriate to support and compute budget. Default
   API values are available choices, not prescribed scientific settings. Small
   groups constrain folds and audit eligibility. Changes to labels or features
   require a fresh evaluation; changing audit/selection settings does not.
5. For an audit, decide whether descriptive contrasts alone answer the question.
   Screening is optional. If requested, explain the chosen mirror ratio and
   eligibility requirement; do not tune them to obtain desired candidate calls.
   For selection, relate delta or population-scale calibration to the user's
   desired granularity. When no scale is established, explore a clearly labeled
   sensitivity path or ask about the target scale rather than inventing certainty.

Use [the interface reference](references/interfaces.md) for CLI configuration,
result reuse and the distinctions needed to interpret scores. Examples there
illustrate syntax; adapt their names and settings to the actual data.

## Execute and communicate

Use `gauge_sc.evaluate`, `audit`, `select`, or their CLI wrappers. Keep input AnnData
unchanged. Save to a new result directory, retaining settings, fingerprints,
validity reasons and tables. Use `gauge.io.load_result` to resume work without refitting.
Use Python when the task needs a custom sklearn classifier or richer Scanpy work.

Show the key table and a figure appropriate to the question. ITCA+ plot helpers
return Matplotlib handles; Scanpy can display explicitly written group labels or
scores on existing embeddings. Export figures only when requested. AnnData
writeback is explicit and requires unchanged cells and original labels.

Explain the outcome in terms of supported distinctions and the specified views:

- Keep ITCA, flat ITCA+ and hierarchical split gain distinct; scores are in nats.
- Positive centered contrast is relative to the eligible-group median. Inspect
  raw contributions and QC covariates before attributing a cause.
- Audit scores are group-level descriptions. Mirror screening is not a p-value,
  FDR guarantee, doublet call, or evidence that a cell should be removed.
- Hierarchy decisions use conditional within-parent scores and ancestor closure;
  different branches can stop at different depths under one global penalty.
- Report ineligible groups, no-threshold outcomes and relevant sensitivity. Treat
  predictability as evidence about a view, not independent biological truth.

Finish with the scientific finding, chosen settings, validity limitations and
paths to useful outputs. Keep operational details proportional to the task.
