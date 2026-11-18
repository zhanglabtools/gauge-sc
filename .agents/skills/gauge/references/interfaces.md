# Interfaces and result reuse

Use the installed version's help and Python signatures when versions differ.
The public computations are `evaluate`, `audit`, and `select`; plotting is under
`gauge.pl`, result loading under `gauge.io`, teaching data under `gauge.datasets`.

## CLI

`gauge_sc` and `python -m gauge_sc` have the same commands:

```bash
gauge inspect data.h5ad
gauge evaluate data.h5ad --config evaluation.json --output results/evaluation
gauge audit results/evaluation --reference molecular --contrast qc --output results/audit
gauge select results/evaluation --view molecular --delta 0.02 0.05 --output results/selection
```

The names, paths and penalties above are syntax examples, not required settings.
`evaluate` reads JSON API keyword arguments. Choose an obs grouping column (or
list of candidates) **or** an ordered hierarchy, plus named views. Each view has
one source:

```json
{
  "groupby": "labels",
  "views": {
    "molecular": {"obsm": "X_pca"},
    "qc": {"obs": ["total_counts", "pct_counts_mt"], "log1p": ["total_counts"]}
  }
}
```

Other sources are `{"X": true}` and `{"layer": "log_normalized"}`.
Optional evaluation fields: `n_splits`, `random_state`, `aggregation`
(`pooled` or `fold_mean`), and `min_child_cells` for hierarchy input.
The CLI uses the API's default classifier; use Python's `classifier=` for a
custom sklearn estimator/pipeline. X, layers and obsm are used as supplied.
Numeric obs views use training-fold standardization and only explicit log1p.

Audit accepts `--groupby` (required with multiple candidates),
`--min-test-cells`, and optional `--maximum-mirror-ratio` in (0, 1]. Omitting
the ratio disables candidate screening. Selection requires exactly one of
`--delta` (nonnegative nats) or `--population-ratio` in (0, 1/e], with one or
more distinct values. Calibration gives delta = -p log(p); it is not a cell-size
filter. Neither command refits a saved evaluation.

Successful analysis commands print a JSON summary and save `result.json` plus
CSV tables. Failures return exit status 2 and a JSON error on stderr. Help and
version output are plain text. Existing output directories are never overwritten.
JSON is authoritative: keep string IDs and boolean flags when loading; do not
reconstruct results from CSV with inferred numeric labels.

## Python and interpretation

```python
import gauge_sc as gauge

evaluation = gauge.io.load_result("results/evaluation")
audit = gauge.audit(evaluation, reference="molecular", contrast="qc")
print(audit.table)
fig = gauge.pl.audit(audit, show=False)
```

Evaluation has `.table` (candidate/view), `.groups` (additive contributions),
`.folds`, `.predictions`, and `.metadata`. Groups use fixed whole-dataset
proportions; views of one grouping share stratified folds. Different groupings
are stratified separately. Fold variation is descriptive, not a confidence interval.

Audit subtracts reference contributions from contrast contributions, then
subtracts the median difference over eligible groups within that grouping.
Review `.metadata['screening']`, eligibility and raw differences. An absent
threshold or no eligible groups is a valid descriptive outcome.

Flat selection retains ITCA, penalty and ITCA+ in `.table`, with selected rows
in `.path`. Hierarchical selection has `.splits` and `.assignments`; top-level
roots remain fixed, ineligible/unary splits block descendants, and the objective
is total conditional split gain rather than flat partition ITCA+.

`.save(new_directory)` exists on every result kind. Loading restores tables,
settings, fingerprints and exclusion reasons; it does not retain fitted models.
`audit.write_to_adata` broadcasts group scores to cells. Hierarchical selection
writeback also needs an explicit chosen penalty when the path has multiple values.
Do not overwrite existing obs/uns keys or treat broadcast scores as cell-level calls.
