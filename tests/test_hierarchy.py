"""Check ancestor closure, conditional normalization and global writeback."""

import itertools
import json

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from sklearn.tree import DecisionTreeClassifier

import gauge_sc as ip


@pytest.fixture
def tree_data():
    leaf = np.repeat(np.arange(8), 20)
    obs = pd.DataFrame(
        {
            "class": np.where(leaf < 4, "A", "B"),
            "sub": [f"s{i//2}" for i in leaf],
            "leaf": [f"l{i}" for i in leaf],
        },
        index=[f"c{i}" for i in range(len(leaf))],
    )
    a = ad.AnnData(leaf[:, None].astype(float), obs=obs)
    a.obsm["spatial"] = np.column_stack([leaf, np.arange(len(leaf)) % 20])
    return a


def evaluate(a, **kwargs):
    return ip.evaluate(
        a,
        hierarchy=["class", "sub", "leaf"],
        views={"molecular": {"X": True}, "spatial": {"obsm": "spatial"}},
        classifier=DecisionTreeClassifier(random_state=0),
        **kwargs,
    )


def test_conditional_scoring_and_parent_folds(tree_data):
    r = evaluate(tree_data, aggregation="fold_mean")
    assert len(r.nodes) == 14 and len(r.table) == 12
    np.testing.assert_allclose(r.table.itca, np.log(2))
    for _, group in r.groups.groupby(["node_id", "view"]):
        assert group.proportion.sum() == pytest.approx(1)
        assert group.contribution.sum() == pytest.approx(np.log(2))
    folds = r.predictions.pivot(
        index=["node_id", "obs_name"], columns="view", values="fold"
    )
    pd.testing.assert_series_equal(folds.molecular, folds.spatial, check_names=False)


def test_dp_matches_exhaustive_search_and_rescues_parent(tree_data, monkeypatch):
    r = evaluate(tree_data)
    # A has a negative immediate gain, but its two children make splitting worthwhile.
    values = {"A": 0.05, "B": 0.05, "s0": 0.69, "s1": 0.69, "s2": 0.01, "s3": 0.01}
    r.table["itca"] = r.table.label.map(values)

    def no_fit(*args, **kwargs):
        raise AssertionError("Pruning must not refit")

    monkeypatch.setattr(DecisionTreeClassifier, "fit", no_fit)
    selected = ip.select(r, view="molecular", population_ratio=0.1)
    score = r.table.query("view == 'molecular'").set_index("node_id")
    parents = r.nodes.set_index("node_id").parent_id
    best = 0
    for bits in itertools.product([False, True], repeat=len(score)):
        accepted = set(score.index[np.array(bits)])
        if any(pd.notna(parents[n]) and parents[n] not in accepted for n in accepted):
            continue
        gain = sum(
            score.loc[n, "itca"]
            - selected.path.delta.iloc[0] * (score.loc[n, "n_children"] - 1)
            for n in accepted
        )
        best = max(best, gain)
    assert selected.path.itca_plus_gain.iloc[0] == pytest.approx(best)
    decisions = selected.splits.set_index("label")
    assert decisions.loc["A", "accepted"] and decisions.loc["A", "immediate_gain"] < 0
    assert decisions.loc["B", "terminal"]
    assert not decisions.loc["s2", "active"]
    assert selected.path.n_groups.iloc[0] == 5
    assert selected.assignments.obs_name.is_unique
    assert set(selected.assignments.obs_name) == set(tree_data.obs_names)


def test_blocked_ancestors_and_zero_ties(tree_data):
    r = evaluate(tree_data)
    r.table["itca"] = -0.1 * np.log(0.1) * (r.table.n_children - 1)
    tied = ip.select(r, view="molecular", population_ratio=0.1)
    assert tied.path.n_groups.iloc[0] == 2
    r.table["itca"] = 0.69
    r.table.loc[r.table.label == "A", ["eligible", "reason"]] = [
        False,
        "insufficient_child_support",
    ]
    blocked = ip.select(r, view="molecular", population_ratio=0.01)
    assert blocked.splits.loc[blocked.splits.label == "A", "terminal"].all()
    assert not blocked.splits.loc[
        blocked.splits.label.isin(["s0", "s1"]), "accepted"
    ].any()


def test_save_reordered_writeback_and_plots(tree_data, tmp_path):
    r = evaluate(tree_data)
    s = ip.select(r, view="molecular", population_ratio=[0.01, 0.3])
    with pytest.raises(ValueError, match="AMBIGUOUS_RATIO"):
        s.write_to_adata(tree_data)
    reordered = tree_data[::-1].copy()
    s.write_to_adata(reordered, population_ratio=0.3)
    expected = (
        s.assignments.query("population_ratio == 0.3").set_index("obs_name").label
    )
    assert (
        reordered.obs.gauge_sc_selection.astype(str).tolist()
        == expected.reindex(reordered.obs_names).tolist()
    )
    reordered.write_h5ad(tmp_path / "tree.h5ad")
    restored = ad.read_h5ad(tmp_path / "tree.h5ad")
    assert (
        json.loads(restored.uns["gauge_sc_selection"]["metadata_json"])["view"]
        == "molecular"
    )
    with pytest.raises(ValueError, match="OUTPUT_EXISTS"):
        s.write_to_adata(reordered, population_ratio=0.3)
    with pytest.raises(ValueError, match="CELL_MISMATCH"):
        s.write_to_adata(tree_data[:-1].copy(), population_ratio=0.3)
    changed = tree_data.copy()
    changed.obs["leaf"] = "changed"
    with pytest.raises(ValueError, match="GROUPING_MISMATCH"):
        s.write_to_adata(changed, population_ratio=0.3)
    for result, name in [(r, "evaluation"), (s, "selection")]:
        result.save(tmp_path / name)
        assert json.loads((tmp_path / name / "result.json").read_text())[
            "kind"
        ].startswith("hierarchy_")
        ip.pl.table(result, show=False)
    ax = ip.pl.hierarchy(s, population_ratio=0.3, show=False)
    assert ax.get_title() == "Selected hierarchy"
    plt.close("all")


def test_invalid_and_ineligible_hierarchy(tree_data):
    bad = tree_data.copy()
    bad.obs.loc["c0", "sub"] = "s3"
    with pytest.raises(ValueError, match="NON_NESTED"):
        evaluate(bad)
    bad = tree_data.copy()
    bad.obs.loc["c0", "leaf"] = None
    with pytest.raises(ValueError, match="MISSING_LABEL"):
        evaluate(bad)
    with pytest.raises(ValueError, match="AMBIGUOUS_GROUPING"):
        ip.evaluate(
            tree_data,
            groupby="class",
            hierarchy=["class", "sub"],
            views={"x": {"X": True}},
        )
    r = evaluate(tree_data, min_child_cells=1000)
    assert not r.table.eligible.any() and r.table.itca.isna().all()
    s = ip.select(r, view="molecular", population_ratio=0.1)
    assert s.path.n_groups.iloc[0] == 2 and s.path.itca_plus_gain.iloc[0] == 0
    unary = tree_data.copy()
    unary.obs["sub"] = unary.obs["class"]
    r = evaluate(unary)
    s = ip.select(r, view="molecular", population_ratio=0.01)
    assert s.path.n_groups.iloc[0] == 2


def test_direct_delta_path_and_view_isolation(tree_data):
    r = evaluate(tree_data)
    s = ip.select(r, view="molecular", delta=[0, 0.2, 1])
    assert s.path.n_groups.tolist() == [8, 8, 2]
    assert s.path.population_ratio.isna().all()
    assert s.metadata["parameterization"] == "delta"
    r.table.loc[r.table.view == "spatial", "itca"] = 999
    again = ip.select(r, view="molecular", delta=[0, 0.2, 1])
    pd.testing.assert_frame_equal(s.path, again.path)
    s.write_to_adata(tree_data, delta=1)
    assert tree_data.obs.gauge_sc_selection.nunique() == 2
    with pytest.raises(ValueError, match="AMBIGUOUS_PENALTY"):
        ip.select(r, view="molecular", delta=0.1, population_ratio=0.1)
    with pytest.raises(ValueError, match="INVALID_DELTA"):
        ip.select(r, view="molecular", delta=-1)
    ip.pl.hierarchy(s, delta=1, show=False)
    ip.pl.table(s, show=False)
    plt.close("all")


def test_deeper_hierarchy_and_branch_counts(tree_data):
    a = tree_data.copy()
    a.obs["deep"] = [f"{leaf}-{i % 2}" for i, leaf in enumerate(a.obs.leaf)]
    r = ip.evaluate(
        a,
        hierarchy=["class", "sub", "leaf", "deep"],
        views={"molecular": {"X": True}},
        classifier=DecisionTreeClassifier(random_state=0),
    )
    s = ip.select(r, view="molecular", delta=[0, 1])
    assert s.path.n_groups.tolist() == [16, 2]
    assert s.assignments.groupby("delta").obs_name.nunique().tolist() == [
        a.n_obs,
        a.n_obs,
    ]
    # A non-binary split pays for two additional terminal groups.
    a = tree_data.copy()
    a.obs["sub"] = a.obs.leaf.map(
        {
            "l0": "s0",
            "l1": "s1",
            "l2": "s2",
            "l3": "s2",
            "l4": "s3",
            "l5": "s3",
            "l6": "s4",
            "l7": "s4",
        }
    )
    r = evaluate(a)
    s = ip.select(r, view="molecular", delta=0.1)
    root = r.table.query("label == 'A' and view == 'molecular'").iloc[0]
    decision = s.splits.query("label == 'A'").iloc[0]
    assert root.n_children == 3
    assert decision.immediate_gain == pytest.approx(root.itca - 0.2)
