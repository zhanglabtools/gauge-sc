# Command-line interface

Use `gauge` after installing the package, or `python -m gauge_sc` in the same
Python environment. `--help` is available on every command. Input data are read
without changing the source file; outputs always require a new directory.

## Inspect inputs

```bash
gauge inspect data.h5ad
```

Returns a JSON inventory of cell/feature counts, observation columns, missing
values, unique values, obsm shapes and layers. This does not fit a model or infer
preprocessing. Inspect provenance and actual values before choosing views.

## Evaluate

```bash
gauge evaluate data.h5ad --config evaluation.json --output results/evaluation
```

`evaluation.json` contains public `evaluate` keyword arguments. For example:

```json
{
  "groupby": "labels",
  "views": {
    "molecular": {"obsm": "X_pca"},
    "qc": {"obs": ["total_counts", "pct_counts_mt"], "log1p": ["total_counts"]}
  },
  "n_splits": 5,
  "random_state": 0,
  "aggregation": "pooled"
}
```

Use a list of grouping columns for candidate comparison. For a hierarchy,
replace `groupby` with `"hierarchy": ["class", "subclass", "supertype"]`;
`min_child_cells` is also configurable. View sources and evaluation parameters
follow the [API guide](api.md). The CLI uses the API's default random forest;
custom sklearn estimators are supported through Python's `classifier=` argument.

## Audit or select without refitting

```bash
gauge audit results/evaluation --reference molecular --contrast qc --output results/audit
gauge select results/evaluation --view molecular --population-ratio 0.01 0.05 --output results/selection
```

Audit operates on a flat evaluation. Use `--groupby` if it contains multiple
candidates. `--min-test-cells` sets eligibility; screening is disabled until
`--maximum-mirror-ratio` is supplied. Its descriptive interpretation is unchanged.

Selection accepts flat or hierarchical evaluations. Supply exactly one of
`--delta` (nonnegative nats per group) or `--population-ratio` in (0, 1/e], with
one or more distinct values. These are example settings, not recommended values
for every dataset. Neither audit nor selection needs the original H5AD file.

## Outputs and errors

Analysis commands save `result.json` and the result's CSV tables, then print a
JSON summary with absolute output paths and row counts. Audit includes screening
status; selection includes its selected path. JSON is authoritative, preserving
string identities, boolean flags and metadata. CSV is for convenient inspection.

```python
import gauge_sc as gauge

evaluation = gauge.io.load_result("results/evaluation")
audit = gauge.audit(evaluation, reference="molecular", contrast="qc")
fig = gauge.pl.audit(audit, show=False)
```

`load_result` accepts a directory or `result.json`, restores any of the five
result kinds, and never trains a model. New saves retain empty table schemas and
dtypes. Older schema-1 saves remain readable with dtype inference from JSON.

Failures return exit status 2 and a JSON error on stderr, for example:

```json
{"error": {"code": "OUTPUT_EXISTS", "message": "choose a new result directory."}}
```

Help and version output are plain text. Figures are displayed or exported through
Python when requested; the CLI does not create them automatically.
