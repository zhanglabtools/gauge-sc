# Python API guide

The root computational API contains `evaluate` (ITCA), `audit` (view-specific ITCA contributions) and `select` (ITCA+). Plotting lives under `gauge_sc.pl`; generated teaching data under `gauge_sc.datasets`; result loading under `gauge_sc.io`.

Use this import for the API examples.

```python
import gauge_sc as gauge
```

## evaluate

```python
evaluate(adata, *, groupby=None, views, n_splits=5, random_state=0,
         aggregation="pooled", classifier=None, hierarchy=None, min_child_cells=10)
```

`groupby` is an obs column or a list of distinct obs columns. Labels may be strings or numbers and are represented as unambiguous strings in result tables. Missing labels and ambiguous string representations are rejected.

Each named view has exactly one source:

```python
{"molecular": {"obsm": "X_pca"}}
{"expression": {"layer": "log_normalized"}}
{"expression": {"X": True}}
{"qc": {"obs": ["total_counts", "pct_counts_mt"], "log1p": ["total_counts"]}}
```

Sources must be finite numeric matrices. The evaluator standardizes obs views within each training fold and uses X and layer views as supplied. Sparse input remains sparse. Use in-memory AnnData and numeric audit covariates.

Custom sklearn classifiers are cloned per fold. Set their random-state parameters explicitly. Only fixed supplied features are evaluated; use a classifier pipeline when additional learned feature processing must occur inside folds.

## audit

```python
audit(result, *, reference, contrast, groupby=None,
      min_test_cells=5, maximum_mirror_ratio=None)
```

One grouping per audit. `groupby` is required if evaluation contains multiple candidates. The raw difference is contrast-view contribution minus reference-view contribution. Its median over eligible groups is subtracted to give `contrast`. The centering scope is exactly this grouping, not a pool of samples, methods or candidate partitions. Independent datasets should be audited separately; paper-specific pooled audit scopes are deferred.

The optional candidate threshold is the smallest observed nonzero absolute centered contrast satisfying `(1 + negative_tail_count) / positive_tail_count <= maximum_mirror_ratio`. Tail counts include equality. A missing threshold means no candidates; disabled screening omits the candidate column entirely. No eligible groups produces an explicit status and missing contrasts, not fabricated zeros.

Results include eligibility reasons, per-view recall and contributions, center, threshold, settings and evaluation metadata. `write_to_adata` requires unchanged labels and the exact original cell set but accepts a different row order.

## Plotting

```python
gauge.pl.audit(
    audit, panel="both", annotate="candidates", sort="contrast",
    figsize=None, palette="vega20", color=None, candidate_color="#B85C70",
    scatter_kws=None, text_kws=None, ax=None, show=True,
    context="paper", font_scale=1, fontfamily=None,
)
gauge.pl.table(result, ax=None, show=True, decimals=3,
            figsize=None, color=None, fontsize=None, palette="vega20",
            context="paper", font_scale=1, fontfamily=None)
```

`pl.audit` returns a Figure. `panel` selects `both`, `contributions` or `contrast`.
Provide one Axes for a single panel or two distinct axes on the same figure for
both. Set the size of your existing figure when passing `ax`; combining `ax`
and `figsize` is an error. The default figure grows with group count.

Scatter labels are `candidates` by default: an unscreened audit has no selected
candidates and therefore no scatter annotations. Use `all`, `none`, or an
explicit list of group names. Labels are shifted to reduce collisions. The
contrast panel always displays group names and orders them by descending
contrast; `sort="input"` preserves the result-table order. Ineligible groups
remain visible without a fabricated zero marker.

For example:

```python
fig = gauge.pl.audit(
    audit,
    annotate=["Q0", "Q1", "Q2"],
    figsize=(9.6, 4.2),
    scatter_kws={"s": 48, "alpha": 0.85},
    text_kws={"fontsize": 9},
    show=False,
)
fig.savefig("audit.svg", bbox_inches="tight")
```

`scatter_kws` applies to points in both panels and accepts `s` (a single positive
area), `alpha`, `marker`, `linewidths`, and `edgecolors`. `text_kws` applies to
scatter annotations and accepts `fontsize`, `fontweight`, `fontstyle`,
`fontfamily`, and `alpha`. Unknown keys and reserved color/coordinate keys raise
an error instead of being ignored. Advanced customizations are possible through
the returned Matplotlib artists. These dictionaries are deliberately scoped,
not unrestricted passthroughs to every artist in a composite figure.

`pl.table` returns an Axes. Scores retain their raw values and integer counts
remain integers. Audit status distinguishes not eligible, not screened,
candidate and not selected. White backgrounds and thin black rules separate rows and headers without
score shading; there is no composite quality ranking.

Plotting uses stored results without modifying them or leaving changes in global Matplotlib
rcParams. Colors identify groups, not the sign of contrast as good/bad. Screening
candidates retain their group colors and receive bold contrast-row labels;
`candidate_color` controls the screening-threshold line. Positive centered contrast is relative to the eligible
group median and need not mean a positive raw difference.

### Context, font size and font family

Both plotting entry points accept `context="paper"`, `font_scale=1`, and
`fontfamily=None`. Presets use `seaborn.plotting_context` directly rather than a
copied approximation. Supported contexts are paper, notebook, talk and poster.
Named contexts keep ITCA+'s black text/axes, white background and outward ticks.
Font family is inherited unless explicitly supplied; availability is system-dependent.

```python
gauge.pl.audit(audit, context="paper", font_scale=1.1, fontfamily="DejaVu Sans")
gauge.pl.table(audit, context="paper", font_scale=1.1, fontfamily="DejaVu Sans")
```

With Seaborn 0.13.2, paper at scale 1 uses 9.6 pt titles/axis labels/body text,
8.8 pt ticks/legend/annotations and 1 pt axes. `font_scale` scales fonts only;
changing the context scales fonts, line widths and marker sizes together.
`text_kws` font sizes and table `fontsize` are explicit point-size overrides.
For example, `text_kws={"fontsize": 8}` stays 8 pt even if font_scale is 1.5.

Set `context=None` to inherit the active Matplotlib/Seaborn context, including
axis colors, spines and grid preferences. `font_scale` can still scale inherited
font sizes, and `fontfamily` can override the inherited family.

```python
with sns.plotting_context("notebook", font_scale=1.2):
    gauge.pl.audit(audit, context=None)
```

An exception-safe local `matplotlib.rc_context` restores settings after every
call. No global theme is installed, and plotting never writes files by default.
Font and context settings also apply to supplied Axes. Result tables share the
same context rather than imposing a separate set of hard-coded point sizes.

### Categorical versus continuous color

`palette="vega20"` uses the 20 color values of Scanpy's
`sc.pl.palettes.vega_20_scanpy` (the adjusted, reordered palette, not the original
`vega_20`). The core package stores the color values without importing Scanpy.
Colors are assigned in result-table input order before plot sorting. For stable
colors across different category orders or subsets, pass a mapping explicitly:

```python
palette = dict(zip(adata.obs["group"].cat.categories,
                   sc.pl.palettes.vega_20_scanpy))
gauge.pl.audit(audit, palette=palette)
gauge.pl.table(audit, palette=palette)
sc.pl.umap(adata, color="group", palette=palette)
```

You can also pass a color sequence, or an existing Scanpy category-color mapping
built from `adata.uns["group_colors"]`. Mapping keys must match result group names
(strings). Missing colors are errors; colors are never silently recycled when
there are more than 20 groups. Supply a longer palette or use `palette=None` for
monochrome. For `pl.audit`, an explicit `color=` overrides the palette.

Unavailable groups remain neutral in audit plots. Vega20 is categorical; a
continuous signed contrast map still needs a diverging colormap such as `RdBu_r`
with its center at zero. Palette choice does not change evaluation or screening.

Default plot and table text is black for readability. Vega20 identifies groups
through data points/lines and small table-row swatches; text colors do not carry
group identity. Axes and outward ticks are black; decorative grids and row stripes are disabled.
Identity and zero reference lines are thin black dashed lines.

### Figure export

The default two-panel audit figure is 7.2 inches wide (about 183 mm), with black
Seaborn paper-context axes (currently 1 pt), outward ticks and compact type. Height grows with group count.
Adjust `figsize` for the intended final placement rather than shrinking a large
figure after export. Defaults are a print-oriented starting point, not a guarantee
of meeting every journal's specifications.

Use a local Matplotlib context for vector export with editable text/embedded fonts:

```python
with plt.rc_context({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"}):
    fig.savefig("audit.pdf", bbox_inches="tight")
    fig.savefig("audit.svg", bbox_inches="tight")
```

This does not alter global plotting settings. SVG text uses the available font
on the viewing system; PDF embeds TrueType fonts.

### Many groups in audit figures

`gauge.pl.audit(audit, groups_per_panel=10)` keeps at most ten rows in each contrast facet. Facets follow one global sort, share contrast limits and row spacing, and wrap at three columns. Every group remains visible; screening and median centering are never recomputed for a facet. Contribution and contrast panels use square plotting areas. Increase `groups_per_panel` to keep more rows together, or use `panel="contrast"` for a dedicated overview. Custom `ax` must supply one axis per contrast facet plus one for contributions if requested.

Above 20 groups, provide a sufficiently large explicit palette (prefer the existing Scanpy group-color mapping), or use `palette=None` for monochrome. Default Vega20 does not silently recycle colors. For very large reports, use a standalone contrast figure and the complete result table rather than shrinking labels to fit a journal-width figure.

## Screening explanations

`audit.metadata["screening"]` adds a machine-readable reason and a plain-language
`message` without changing the screening rule. The existing `screening_status`
remains `disabled`, `no_threshold` or `threshold_found` for compatibility.

| Reason | Meaning |
|---|---|
| `disabled` | No ratio requested; continuous contrasts only. |
| `no_eligible_groups` | No groups meet held-out support requirements. |
| `insufficient_positive_groups` | Too few strictly positive eligible contrasts, even assuming no negative-tail groups. |
| `mirror_ratio_not_met` | Enough positive groups exist, but no observed threshold satisfies the tail-ratio rule. |
| `threshold_found` | A qualifying threshold and candidate set were found. |

Additional fields are `eligible_count`, `positive_group_count`,
`minimum_positive_count`, `threshold`, `positive_tail_count`,
`negative_tail_count`, `achieved_mirror_ratio` and `candidate_count`.
Positive-group count counts all strictly positive eligible contrasts; tail counts
refer only to the chosen threshold, with equality included. Without a threshold,
tail counts and achieved ratio are `None`; candidate count is zero. When screening
is disabled, the minimum required count is also `None` and the result table has no
candidate column. These explanations survive JSON exports and AnnData writeback.

Default screening remains disabled (`maximum_mirror_ratio=None`). For explicit
exploratory screening we recommend starting at `1 / 3`, comparing 0.2 and 0.5,
and reporting all checked settings. This is an operational recommendation, not a
calibrated error rate or a guarantee of suitable power. `1 / 3` needs at least
three positive-tail groups; `0.33` needs four. Do not round the supplied parameter.

## select

```python
select(result, *, view, population_ratio=None, delta=None)
```

Pass one `EvaluationResult` and explicitly name the selection view (normally molecular). `population_ratio` is a positive scalar or a list of distinct values up to `1/e`. There is no implicit default. This monotone interval is an API constraint: larger ratios mean larger penalties. With `p = population_ratio`, `delta = -p * log(p)` and `penalty = delta * n_groups`, the criterion is `itca_plus = itca - penalty`. Natural logs give scores in nats. This matches the research implementation formerly named ITCA-MDL; the release name is ITCA+.

Population ratio calibrates the complexity penalty, not minimum observed group size. No groups are deleted, no classifier is retrained, and no new candidate is generated. Existing evaluation aggregation (pooled or fold mean) is retained and recorded. Use `aggregation="fold_mean"` when that is the intended manuscript convention; selection does not silently change it. Candidates must have complete predictions on the same observations in the same evaluation/view. Different candidates retain their own stratified folds; selection does not turn them into paired-fold estimates.

`SelectionResult.table` contains every candidate at every requested ratio, with explicit `itca`, `delta`, `penalty`, `itca_plus` and `selected` columns. `path` contains one selected row per ratio in the requested order. Highest ITCA+ wins; exact ties prefer fewer groups, then lexical grouping-column order. No near-tie tolerance or confidence claim is implied. `save(new_directory)` exports both tables and evaluation provenance; `gauge.pl.table(selection)` displays the three distinct quantities. For long paths, inspect `selection.path` rather than putting all rows into a single publication figure.

For flat EvaluationResult input, selection compares supplied candidates. HierarchyEvaluationResult input instead uses the conditional split objective described below. ITCA+ is a candidate-level objective; group-level contributions returned by `evaluate` and used by `audit` remain ITCA contributions.

## Hierarchy evaluation and selection

```python
evaluation = gauge.evaluate(
    adata, hierarchy=["class", "subclass", "supertype"],
    views={"molecular": {"X": True}, "spatial": {"obsm": "spatial"}},
    aggregation="fold_mean", min_child_cells=10,
)
selection = gauge.select(evaluation, view="molecular", delta=[0.1, 0.4, 0.8])
selection.write_to_adata(adata, delta=0.4)
gauge.pl.hierarchy(selection, delta=0.4)
```

`hierarchy` and `groupby` are mutually exclusive. Supply at least two distinct
nested obs columns, ordered coarse to fine, with complete labels. Within a level,
each child label must have exactly one parent; non-nested labels are rejected.
Nodes have unambiguous JSON path identifiers. The first level forms a fixed forest
of roots; the algorithm does not score merging these roots or a synthetic root.

For each internal node, prediction occurs only among its cells, using its immediate
child labels. Child proportions are normalized within that parent. Views share
folds within the parent. Representations remain fixed; no local PCA is silently
refitted. At least two children and at least `max(min_child_cells, n_splits)` cells
per child are required. Otherwise scores are missing with an eligibility reason.
Ineligible splits, including unary nodes, block descendants, matching the research
rule; nodes/cells are not discarded. Reorganize a taxonomy explicitly if unary
levels should instead be bypassed.

`HierarchyEvaluationResult` retains `table`, `nodes`, `assignments`, `groups`,
`folds`, `predictions` and per-node evaluation provenance. `save` exports these.
It is not input for mirror-tail `audit()`: inspect conditional view scores for the
molecular-selected splits separately, rather than pooling unlike parent populations.

For node v with k children, the split gain is `conditional_itca(v) - delta*(k-1)`.
Dynamic programming compares stopping (gain zero) with splitting plus the best
allowed gains of its children. A negative immediate split can be worthwhile if
its descendants compensate. Ineligible ancestors block all descendants. Exact
zero subtree gain favors stopping. This is globally optimal for this fixed tree
and the stored additive conditional scores; it is not a claim about the biological
truth or arbitrary groupings outside that tree.

The selected path reports `conditional_itca_sum`, `penalty`, `itca_plus_gain`,
`n_groups` (terminals) and `n_splits` (accepted internal nodes, not CV folds).
Penalty is `delta * (n_terminals - n_roots)`: these are gains relative to the fixed
roots. **The conditional sum is not a re-evaluation of flat terminal labels.**
To obtain their flat ITCA, explicitly call `evaluate(groupby="gauge_sc_selection", ...)`
after writeback; it need not equal the conditional sum.

`HierarchySelectionResult.splits` contains every node and parameter with active,
terminal, accepted, eligible, reason, immediate_gain and subtree_gain fields. A
positive subtree gain at an unreachable descendant is only a hypothetical local
optimum; `accepted` says whether it is selected globally. `assignments` covers each
evaluated observation exactly once per parameter. Spatial/QC views never influence
selection unless explicitly named as `view`.

Writeback requires the full unchanged cell set and hierarchy, allows row reordering,
rejects AnnData views and existing output keys, and stores a categorical terminal
label plus serializable metadata. Select a single stored `delta` or `population_ratio`
for a multi-point path. After writeback, a spatial section can be displayed as a
subset without recomputing its pruning. `pl.hierarchy` returns Axes, uses the same
local context/font controls as other plots, and never saves automatically. Colors
can be supplied as a mapping from `level: label` to match the Scanpy spatial map.

### Direct penalties

Both flat and hierarchy selection accept a scalar or distinct sequence `delta >= 0`
as an alternative to `population_ratio`. Supply exactly one parameterization.
With population ratios, `delta=-p log(p)`; `p` remains restricted to `(0,1/e]` so
increasing it increases the cost. Direct delta is needed for paths above `1/e`,
including stronger penalties in the research hierarchy analyses. No population
ratio is invented for those values: its output field is null. `delta=0` disables
the complexity penalty. Neither parameter filters small groups. In hierarchy
mode, population-scale interpretation is conditional on each parent, not a global
minimum cell fraction.

## Portable table typography

Use `gauge.pl.table(result, fontfamily="Inter")` for the bundled Inter typeface; registration is local to the Python process and needs no font download or OS installation. The default `fontfamily=None` still respects the user's Matplotlib environment. Font licensing is included with the package. In notebooks and webpages, prefer `display(result.table.round(3))` for selectable, accessible HTML tables; exporting a figure remains optional.

## Real tutorial datasets

`gauge.datasets.load_example(name, data_dir=None, cache_dir=None, download=True)` returns one verified AnnData. See [dataset names, public downloads and offline setup](real-data-tutorials.md). The computational API contains `evaluate`, `audit` and `select`.

## Saving and resuming

Every returned result has `save(new_directory)`, which writes versioned JSON and
CSV tables without overwriting. Use `gauge.io.load_result(directory_or_result_json)`
to restore tables, metadata, identities and validity reasons. JSON is authoritative;
new saves also retain dtypes and empty table schemas. No fitted estimator is stored.

```python
evaluation.save("results/evaluation")
restored = gauge.io.load_result("results/evaluation")
audit = gauge.audit(restored, reference="molecular", contrast="qc")
selection = gauge.select(restored, view="molecular", delta=[0.02, 0.05])
```

Audit and selection reuse stored evaluation evidence. Labels, input views and
evaluation settings must be reevaluated if changed. The [CLI](cli.md) exposes the
same computations; [Work with your agent](agents.md) explains the installable skill.
