# Tutorial data policy

The Python distribution contains no real dataset and downloads nothing on installation or import. The first notebook uses deterministic, explicitly synthetic teaching data generated in memory. This exercises the user workflow but is not a validation of manuscript findings.

Real-data tutorials use small, fixed, documented subsets from the [GAUGE data release](https://github.com/zhanglabtools/gauge-sc/releases/tag/tutorial-data-v1). The release documents each dataset's source, redistribution terms, processing, sampling and dimensions. `gauge.datasets.load_example()` automatically verifies the pinned SHA-256 of each file it loads. See the [tutorial data guide](real-data-tutorials.md) for working with downloaded, cached or local files. Dataset files are distributed separately from Git checkouts and Python packages.

Keep original group proportions when possible: equalizing groups changes ITCA information weights. Explain any enrichment of rare groups. Recompute scores after subsampling. Keep QC provenance explicit when a full gene set is reduced. Do not choose subsets or seeds by searching for a desired audit result.

Testing uses generated miniature inputs, so core tests require no network or manuscript datasets. Tutorial notebooks contain code, prose and compact output, never embedded h5ad data. Real application notebooks retain executed outputs and a source hash; synthetic notebooks are executed during each documentation build.

## Controlled audit example

`audit_example()` creates seven equal-sized groups (80 cells each; 560 cells, 50 features). M0–M3 have distinct molecular marker blocks and share a QC distribution. Q0–Q2 share a fifth molecular program but have distinct, independently simulated QC distributions. QC-like column names support the standard workflow, but these covariates are **not computed from X**. X contains nonnegative synthetic features, not counts; there is no counts layer. `simulation_role` is explanatory metadata, never a predictor.

The tutorial enables mirror screening at `maximum_mirror_ratio=1/3`. The +1 numerator requires at least three positive-tail groups even with no negative-tail groups; four molecular groups provide a majority for median centering. Equal sizes are not a recommendation to rebalance real data.

Seeds 0–4 form a fixed regression check, not a benchmark or evidence of error-rate control. Positive centered contrast, a screening candidate, and known simulation roles are different concepts. Real data may legitimately yield no candidates.

## Additional application tutorials

The grouping-comparison notebook fixes three five-group assignments on the same
560 synthetic observations: a simulated molecular reference, expression KMeans
and QC KMeans. QC covariates keep the independent simulation provenance described
above; clustering is performed once before evaluation.

The transferred-label notebook generates source and target with fixed independent
seeds 0 and 1, gives them distinct observation identifiers, and transfers the five
source molecular labels using nearest neighbors in the common synthetic feature
space. The source and target are deliberately compatible teaching samples, not a
simulation of real cross-platform domain shift. Target simulation assignments
are removed before evaluation. The evaluator receives target observations only;
no real dataset or fitted transfer model is shipped with the package.

## Stress-check inputs

The controlled stress example reuses the five molecular programs from the audit
simulation. Label permutations preserve the complete label histogram. Cell
subsets are nested prefixes of a fixed uniformly shuffled observation order;
no label balancing or post hoc subset selection is used. Expression/PCA is kept
fixed to isolate conditional scoring behavior. Seeds 0–2 act on a fixed perturbed
label vector. All modified AnnData objects are copies; no files are downloaded
or automatically exported.

## Hierarchy and toy spatial example

The seventh notebook builds a fixed three-level taxonomy from the existing seven
synthetic groups: two root classes, four subclasses and seven leaf labels. One
subclass has a single child and is explicitly ineligible for splitting. Toy 2D
coordinates are sampled around predefined label-associated centers with seed 0;
this is neither realistic tissue simulation nor evidence of anatomical validity.
All 560 observations are evaluated globally. Molecular and spatial views are
scored separately, and only the molecular view enters pruning. Color-matched tree
and map display the same stored global terminal assignments.
