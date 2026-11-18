# Numerical reproducibility

Small differences in preprocessing or held-out predictions can change the
selected grouping when candidate scores are close. Report the ITCA+ margin
between candidates to show how sensitive the choice is.

The [Allen tutorial](../notebooks/real_data/14_allen_granularity.ipynb) displays the two candidates' raw ITCA, total complexity cost and penalized ITCA+ at `δ=0.02`, with their difference reported to nine decimal places.
Do not round scores before selection, change the penalty to recover a preferred
answer, or describe a fixed random seed as cross-platform bitwise reproducibility.

Run the independent fresh-process audit from an installed tutorial environment:

```bash
python scripts/check_allen_stability.py --data /path/to/tutorial-data \
  --output artifacts/allen-stability --coretypes default Haswell --seeds 0 1 2
```

Use `Haswell` only on a compatible x86/OpenBLAS environment. The audit records the actual architecture,
library versions and threads reported by threadpoolctl; inspect these to confirm
whether a requested core type took effect. Seeds vary the evaluation folds and
classifier; PCA and candidate generation remain fixed to isolate that comparison.
For another platform use `--coretypes default`.

Each run preserves PCA, labels, input fingerprints, folds, predictions, full
precision scores and margins. Compare partition assignments and biological
readouts as well as group count. An environment lock provides a reference run;
portability still requires explicit numerical and decision sensitivity checks.
This audit characterizes sensitivity, and does not establish biological truth.
