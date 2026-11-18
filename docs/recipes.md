# Optional recipes

Keep the introductory notebook focused on the concept. Details and optional
steps live here and in [API reference](api.md).

Import GAUGE before running the recipes.

```python
import gauge_sc as gauge
```

## Candidate screening

```python
screened = gauge.audit(result, reference="molecular", contrast="qc",
                    maximum_mirror_ratio=1/3)
print(screened.metadata["screening"]["message"])
print(screened.table)
```

This applies a descriptive mirror-tail rule, not formal FDR control. No qualifying
threshold means no candidates, which is a valid outcome. It does not refit the
classifier. The default continuous audit omits candidate calls entirely.

## Screening sensitivity without retraining

Use `None` for the default continuous audit. For exploratory candidate screening,
start at `1 / 3` and compare prespecified settings; it is not an FDR target.
The +1 numerator requires at least `ceil(1/q)` positive-tail groups even when
there are no negative-tail groups. `0.33` is stricter than `1 / 3` and requires
four rather than three positive-tail groups.

```python
import pandas as pd

rows = []
for ratio in [0.2, 1 / 3, 0.5]:
    audit = gauge.audit(result, reference="molecular", contrast="qc",
                     maximum_mirror_ratio=ratio)
    info = audit.metadata["screening"]
    rows.append({
        "maximum_mirror_ratio": ratio,
        "threshold": info["threshold"],
        "positive_tail": info["positive_tail_count"],
        "negative_tail": info["negative_tail_count"],
        "achieved_ratio": info["achieved_mirror_ratio"],
        "candidates": audit.table.loc[audit.table.candidate, "group"].tolist(),
        "reason": info["reason"],
    })
comparison = pd.DataFrame(rows)
```

Every call uses the same EvaluationResult; no model is refitted. In the seven-group
example, 0.2 has too few positive groups; 1/3 and 0.5 select Q0–Q2. This agreement
is parameter stability, not evidence that candidates are technical artifacts.
No threshold means tail counts are undefined (`None`), not zero.
See the [executable sensitivity notebook](../notebooks/02_audit_sensitivity.ipynb).

## Styling

For uniform font and line scaling, use `context="paper"`, `font_scale=1.1` and
optionally `fontfamily="DejaVu Sans"` on either plotting function. Use
`context=None` inside your existing Seaborn context to inherit its settings.


```python
fig = gauge.pl.audit(audit, panel="contrast", figsize=(6, 4),
                  scatter_kws={"s": 52, "alpha": 0.85}, show=False)
fig.axes[0].set_title("Relative QC contribution", loc="left")
fig.savefig("audit.svg", bbox_inches="tight")
```

Use `annotate="all"`, `"none"`, `"candidates"` or explicit group names for scatter
labels. Scoped styling parameters and axis handling are listed in the API guide.
Do not use the categorical Vega20 palette to encode continuous contrast values.

## Export

```python
result.save("outputs/evaluation")
audit.save("outputs/audit")
adata.write_h5ad("outputs/audited.h5ad")
```

Result directories must not already exist. JSON and CSV exports include actual
observation identifiers, labels and evaluation provenance; they are local outputs,
not files to commit automatically. Explicit AnnData writeback must occur before
saving if broadcast scores are wanted. Exporting and replotting do not retrain.

## Evaluation assumptions

Candidate assignments and supplied representations are fixed before evaluation.
This is held-out prediction of those assignments, not validation of the entire
clustering workflow. Views share stratified folds within one grouping. Numeric
QC preprocessing is fitted inside training folds; supplied PCA is not refitted.

Default recall is pooled over OOF predictions; `aggregation="fold_mean"` is an
explicit alternative. Information weights use fixed observed group proportions.
Audit centers differences over eligible groups within the selected grouping;
positive centered contrast does not necessarily mean a positive raw difference.

Evaluation needs at least `n_splits` observations per group. Audit additionally
requires `min_test_cells` observations in every held-out fold. Ineligible groups
remain visible with reasons, without fabricated zero contrasts.

## Real-data subsets

Subsampling changes group proportions and information weights. Record the policy,
and recompute scores for the sampled cells.
QC covariates should retain their intended full-expression provenance when only a
subset of genes is distributed. See [data policy](data.md).

## Comparing existing groupings

Pass several `.obs` columns to `evaluate(groupby=[...])` with explicitly named
views. Evaluation returns one ITCA score per grouping/view and uses the same
observations for every candidate. Views share folds within a grouping; distinct
groupings have separately stratified folds. Candidate clustering is external and
fixed before evaluation. Prediction support does not independently validate the
candidate-generation procedure or biological annotation. ITCA information weights
depend on group proportions; a larger score alone is not a complexity-calibrated
choice. Use `select` for ITCA+ comparison on one explicit view.

## Transferred labels

```python
# target already contains transferred labels and prepared target features
result = gauge.evaluate(
    target, groupby="transferred_labels",
    views={"molecular": {"obsm": "X_pca"}},
    aggregation="fold_mean",
)
```

Label transfer is upstream. ITCA evaluates fixed transferred assignments using
cross-validation **within the target**, with target group proportions. It does
not reuse the source-fitted classifier as the evaluation model and is not
source-to-target accuracy. cITCA identifies this assignment provenance, not a
separate score formula or API. The tutorial's neighbor classifier is only a
small external demonstration; no transfer algorithm is added to ITCA+.

Real transfer workflows must align features and preprocessing, decide how to
handle unknown/missing labels and explicitly define any shared-label subset.
The library rejects missing labels instead of silently excluding cells; it does
not decide that an unknown label is a biological population. Source and target
scores with different populations should not be combined as though they were
comparable candidate rows in one `select` call. If labels are derived from the
same target features used for scoring, high predictability is still conditional
on that fixed workflow, not independent validation of label correctness.

## Controlled stress checks

The [advanced stress notebook](../notebooks/05_stress_checks.ipynb) changes one
factor at a time. Its label-shuffle fractions 0, 0.25 and 0.5 are prespecified;
permuting within the subset preserves group counts, so the entropy ceiling is
constant. The actual fraction of changed assignments can be smaller because a
permutation may assign the original label. No result-dependent seed search is
performed.

Nested uniform samples of 560, 140 and 20 cells use a fixed order. Subsampling
recomputes group proportions, and the library checks per-group CV support.
Unavailable evaluation is recorded explicitly rather than reported as zero or
repaired by silently dropping/merging labels. A group with five total observations
can support fivefold evaluation but still fail audit's default five-per-test-fold
requirement. The notebook reports both outcomes separately. PCA remains fixed;
this is not a stress test of the whole expression-processing pipeline.

Seeds 0, 1 and 2 use exactly the same 25%-shuffle assignments and representation;
both CV splits and default forest randomness vary. The three resulting scores
are descriptive sensitivity, not confidence limits. Classifier sensitivity is a
separate optional check: pass cloneable estimators via `classifier`, hold
`random_state` fixed for matched splits, set each estimator's own seed explicitly,
and fit any learned transformations in a classifier pipeline. Do not pool
classifier, seed and input-perturbation changes as though they were the same
experiment.

ITCA scores from different cell subsets are not directly reusable as candidates
in a single ITCA+ selection result. The observed population and information
weights have changed. Higher predictability does not establish biological truth.

## Hierarchy and spatial selection

Use `hierarchy=[...]` for nested annotations and explicitly select the molecular
view. `select(evaluation, population_ratio=...)` uses a population-scale penalty;
`select(evaluation, delta=...)` directly sets per-additional-group cost in nats.
Do not provide both. The manuscript's direct-delta paths can exceed the maximum
cost attainable by the monotone population-ratio calibration.

Spatial coordinates can be evaluated as a separate numeric view and displayed
with Scanpy after explicit global-pruning writeback. The built-in example uses
toy coordinates, not a tissue image or anatomical validation. The molecular
objective never automatically incorporates spatial/QC predictability.

To inspect control views for the retained splits, without changing selection:

```python
accepted = selection.splits.query("delta == 0.4 and accepted").node_id
controls = evaluation.table[evaluation.table.node_id.isin(accepted)]
controls.pivot(index="node_id", columns="view", values="itca")
```

These scores are conditional on different parent populations. Do not interpret
their sum as flat ITCA, or run a pooled mirror test on them. Categorical control
features still require explicit upstream encoding; imputation is not automatic.
