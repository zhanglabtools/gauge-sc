"""Optional developer check against a local, separate research repository.

No legacy import or dependency is introduced into the installed package. This
check extracts only the old mathematical definitions for numerical comparison.
It checks the historical pooled score and archived mirror thresholds, not the
full manuscript analysis or identical fitted models.
"""

import argparse
import ast
from pathlib import Path

import numpy as np
import pandas as pd

import gauge_sc as ip
from gauge_sc._audit import _mirror_threshold

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("research_repo", type=Path)
args = parser.parse_args()
root = args.research_repo
namespace = {"np": np}
for path, names in [
    (root / "src/scit/utils.py", {"compute_y_dist", "bidict"}),
    (root / "src/scit/metrics/base.py", {"itca"}),
]:
    tree = ast.parse(path.read_text())
    selected = [
        node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name in names
    ]
    exec(
        compile(ast.Module(body=selected, type_ignores=[]), str(path), "exec"),
        namespace,
    )
adata = ip.datasets.audit_example()
result = ip.evaluate(
    adata,
    groupby="group",
    views={"molecular": {"obsm": "X_pca"}, "qc": {"obs": ["total_counts"]}},
)
errors = []
for row in result.table.itertuples():
    predictions = result.predictions[result.predictions.view == row.view]
    categories = sorted(predictions.group.unique())
    lookup = {label: i for i, label in enumerate(categories)}
    y = predictions.group.map(lookup).to_numpy()
    pred = predictions.predicted_group.map(lookup).to_numpy()
    mapping = namespace["bidict"]({i: i for i in range(len(categories))})
    legacy = namespace["itca"](y, pred, mapping)
    errors.append(abs(legacy - row.itca))
    np.testing.assert_allclose(legacy, row.itca, atol=1e-12, rtol=0)
print(
    f"Pooled ITCA agrees with legacy formula for {len(errors)} views; maximum absolute error {max(errors):.3g}."
)
path = (
    root
    / "results/ding2020/contrast_score_attribution_audit_2026-07-31/all_cluster_mirror_audit.csv"
)
archived = pd.read_csv(path)
count = 0
for panel, frame in archived.groupby("panel_id"):
    values = frame.centered_attribution_contrast.to_numpy()
    threshold = _mirror_threshold(values[np.isfinite(values)], 1 / 3)
    saved = frame.mirror_threshold.dropna().unique()
    if len(saved):
        assert threshold is not None
        np.testing.assert_allclose(saved, threshold, rtol=0, atol=1e-12)
    else:
        assert threshold is None
    actual = (
        np.zeros(len(frame), dtype=bool) if threshold is None else values >= threshold
    )
    np.testing.assert_array_equal(actual, frame.mirror_candidate.to_numpy())
    count += 1
print(
    f"Mirror thresholds and candidate calls agree for {count} archived panels ({len(archived)} groups)."
)
print(
    "Archived panels pool eligible groups across methods. This checks the threshold primitive; the public audit API deliberately centers within one grouping."
)

# ITCA+ is the renamed flat ITCA-MDL criterion, with no legacy runtime import.
from typing import Iterable

source = root / "scripts/exp03_resolution_selection/pipeline.py"
tree = ast.parse(source.read_text())
nodes = [
    node
    for node in tree.body
    if isinstance(node, ast.FunctionDef)
    and node.name in {"compute_mdl_path", "compute_population_scale_table"}
]
selection_namespace = {"pd": pd, "np": np, "Iterable": Iterable}
exec(
    compile(ast.Module(body=nodes, type_ignores=[]), str(source), "exec"),
    selection_namespace,
)
group = adata.obs["group"].astype(str)
adata.obs["merged"] = group.where(~group.str.startswith("Q"), "Q")
evaluation = ip.evaluate(
    adata,
    groupby=["group", "merged"],
    views={"molecular": {"obsm": "X_pca"}},
    aggregation="fold_mean",
)
ratios = [0.001, 0.01, 0.05, 0.1, 0.3]
selection = ip.select(evaluation, view="molecular", population_ratio=ratios)
calibration = selection_namespace["compute_population_scale_table"](adata.n_obs, ratios)
legacy = evaluation.table.rename(columns={"itca": "aITCA", "n_groups": "K"}).assign(
    track="leiden"
)
legacy_path = selection_namespace["compute_mdl_path"](legacy, calibration.delta)
legacy_path = legacy_path.loc[legacy_path.selection_group == "leiden"]
np.testing.assert_allclose(
    legacy_path.mdl_score, selection.path.itca_plus, rtol=0, atol=1e-12
)
assert legacy_path["groupby"].tolist() == selection.path["groupby"].tolist()
print(
    "ITCA+ formula and selected path match research ITCA-MDL at five population ratios."
)

# Ancestor-constrained DP parity with the original three-level research routine.
source = root / "scripts/exp04_spatial_transfer_citca_mdl/pipeline.py"
nodes = [
    node
    for node in ast.parse(source.read_text()).body
    if isinstance(node, ast.FunctionDef)
    and node.name in {"_split_gain", "compute_adaptive_pruning_path_dp"}
]
import json

hierarchy_namespace = {"pd": pd, "np": np, "Iterable": Iterable, "json": json}
exec(
    compile(ast.Module(body=nodes, type_ignores=[]), str(source), "exec"),
    hierarchy_namespace,
)
adata.obs["class"] = np.where(group.str.startswith("Q"), "Q", "M")
adata.obs["subclass"] = group.map(
    {
        "M0": "M01",
        "M1": "M01",
        "M2": "M23",
        "M3": "M23",
        "Q0": "Q0only",
        "Q1": "Q12",
        "Q2": "Q12",
    }
)
adata.obs["leaf"] = group
hierarchy_evaluation = ip.evaluate(
    adata,
    hierarchy=["class", "subclass", "leaf"],
    views={"molecular": {"obsm": "X_pca"}},
    aggregation="fold_mean",
)
records = []
for row in hierarchy_evaluation.table.itertuples():
    children = hierarchy_evaluation.nodes.loc[
        hierarchy_evaluation.nodes.parent_id == row.node_id, "label"
    ]
    records.append(
        dict(
            parent_level=row.level,
            child_level="subclass" if row.level == "class" else "leaf",
            parent_label=row.label,
            eligible=row.eligible,
            split_ITCA_mean=row.itca,
            K_children=row.n_children,
            delta_K=row.n_children - 1,
            child_mapping=json.dumps({label: i for i, label in enumerate(children)}),
        )
    )
deltas = [0, 0.2, 0.4, 0.8]
legacy_path, legacy_splits = hierarchy_namespace["compute_adaptive_pruning_path_dp"](
    pd.DataFrame(records), deltas, adata.n_obs
)
new = ip.select(hierarchy_evaluation, view="molecular", delta=deltas)
np.testing.assert_allclose(
    legacy_path.mdl_objective, new.path.itca_plus_gain, rtol=0, atol=1e-12
)
np.testing.assert_array_equal(legacy_path.total_leaves, new.path.n_groups)
for delta in deltas:
    old = legacy_splits.loc[legacy_splits.delta == delta]
    current = new.splits.loc[(new.splits.delta == delta) & new.splits.accepted]
    assert set(zip(old.parent_level, old.parent_label)) == set(
        zip(current.level, current.label)
    )
print(
    "Hierarchy objectives, terminal counts and accepted splits match research DP at four deltas."
)
