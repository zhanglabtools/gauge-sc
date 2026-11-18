"""Scientific scope checks for the runnable application examples."""

import json
from pathlib import Path

import pandas as pd

import gauge_sc as ip

ROOT = Path(__file__).resolve().parents[1]


def execute_example(name, monkeypatch):
    """Execute actual notebook cells; full rendering is checked separately."""
    monkeypatch.setattr(ip.pl, "table", lambda *args, **kwargs: None)
    notebook = json.loads((ROOT / "notebooks" / name).read_text())
    namespace = {}
    for cell in notebook["cells"]:
        if cell["cell_type"] == "code":
            exec(compile("".join(cell["source"]), name, "exec"), namespace)
    return namespace


def test_comparison_tutorial_uses_same_cells_and_explicit_views(monkeypatch):
    ns = execute_example("03_compare_groupings.ipynb", monkeypatch)
    adata, result = ns["adata"], ns["result"]
    assert len(result.table) == 6
    assert set(result.table.n_groups) == {5}
    assert set(result.table.n_cells) == {adata.n_obs}
    for _, rows in result.predictions.groupby(["groupby", "view"]):
        assert set(rows.obs_name) == set(adata.obs_names)
        assert rows.obs_name.is_unique
    scores = result.table.pivot(index="groupby", columns="view", values="itca")
    assert scores.loc["reference", "molecular"] > scores.loc["qc_clusters", "molecular"]
    assert scores.loc["qc_clusters", "qc"] > scores.loc["reference", "qc"]


def test_transfer_tutorial_evaluates_target_not_source(monkeypatch):
    ns = execute_example("04_transferred_labels.ipynb", monkeypatch)
    source, target, result = ns["source"], ns["target"], ns["result"]
    assert not set(source.obs_names) & set(target.obs_names)
    assert set(result.predictions.obs_name) == set(target.obs_names)
    assert not {"group", "simulation_role"} & set(target.obs.columns)
    assert set(result.table["groupby"]) == {"transferred"}
    assert result.metadata["classifier_class"].endswith("RandomForestClassifier")
    assert result.metadata["aggregation"] == "fold_mean"
    groups = result.groups.query("view == 'molecular'").set_index("group")
    counts = target.obs.transferred.value_counts().reindex(groups.index)
    assert groups.n_cells.tolist() == counts.tolist()
    assert groups.proportion.tolist() == (counts / target.n_obs).tolist()
    folds = result.predictions.pivot(index="obs_name", columns="view", values="fold")
    pd.testing.assert_series_equal(folds.molecular, folds.qc, check_names=False)


def test_stress_tutorial_preserves_controlled_factors_and_reports_failure(monkeypatch):
    import matplotlib.pyplot as plt
    import numpy as np

    monkeypatch.setattr(plt, "show", lambda: None)
    ns = execute_example("05_stress_checks.ipynb", monkeypatch)
    adata = ns["adata"]
    expected = adata.obs.labels.value_counts().sort_index()
    try:
        for sample in ns["label_inputs"].values():
            pd.testing.assert_series_equal(
                sample.obs.labels.value_counts().sort_index(), expected
            )
            np.testing.assert_array_equal(sample.obsm["X_pca"], adata.obsm["X_pca"])
        np.testing.assert_array_equal(adata.obs.labels, ns["original_labels"])
        label_summary = ns["label_summary"]
        assert (
            label_summary["Actually changed"] <= label_summary["Shuffled subset"]
        ).all()
        np.testing.assert_allclose(label_summary.Entropy, label_summary.Entropy.iloc[0])
        assert label_summary.ITCA.iloc[-1] < label_summary.ITCA.iloc[0]
        sample_summary = ns["sample_summary"].set_index("Cells")
        assert pd.isna(sample_summary.loc[20, "ITCA"])
        assert sample_summary.loc[20, "Status"] == "Insufficient cells per group for CV"
        assert 20 not in ns["sample_results"]
        assert (
            sample_summary.loc[140, "Audit-eligible groups"]
            < sample_summary.loc[140, "Groups present"]
        )
        assert set(ns["sample_inputs"][20].obs_names) <= set(
            ns["sample_inputs"][140].obs_names
        )
        for count, result in ns["sample_results"].items():
            sample = ns["sample_inputs"][count]
            assert set(result.predictions.obs_name) == set(sample.obs_names)
            expected_p = sample.obs.labels.value_counts(normalize=True)
            observed = result.groups.set_index("group").proportion
            np.testing.assert_allclose(observed, expected_p.reindex(observed.index))
        first = ns["seed_results"][0]
        for seed, result in ns["seed_results"].items():
            assert result.metadata["random_state"] == seed
            assert (
                result.metadata["grouping_fingerprints"]
                == first.metadata["grouping_fingerprints"]
            )
            assert (
                result.metadata["view_fingerprints"]
                == first.metadata["view_fingerprints"]
            )
    finally:
        plt.close("all")


def test_concept_tour_demonstrates_information_cost_and_view_dependence(monkeypatch):
    import matplotlib.pyplot as plt
    import numpy as np

    monkeypatch.setattr(plt, "show", lambda: None)
    ns = execute_example("00_concepts.ipynb", monkeypatch)
    try:
        table = ns["comparison"]
        assert table.loc["coarse", "held_out_accuracy"] == 1
        assert table.loc["resolved", "held_out_accuracy"] == 1
        assert table.loc["resolved", "itca"] > table.loc["coarse", "itca"]
        assert ns["selection"].path["groupby"].tolist() == ["resolved", "coarse"]
        guide = (ROOT / "docs/concepts.md").read_text()
        for name, row in table.iterrows():
            snapshot = f"| {name.title()} | {int(row.n_groups)} | {row.held_out_accuracy:.3f} | {row.itca:.3f} |"
            assert snapshot in guide
        for row in ns["selection"].path.itertuples():
            snapshot = f"| {row.population_ratio:.2f} | {row.groupby.title()} · {row.n_groups} groups | {row.itca:.3f} | {row.penalty:.3f} | {row.itca_plus:.3f} |"
            assert snapshot in guide
        for row in ns["result"].table.itertuples():
            groups = ns["result"].groups.query(
                "groupby == @row.groupby and view == @row.view"
            )
            np.testing.assert_allclose(
                groups.contribution,
                -groups.proportion * np.log(groups.proportion) * groups.recall,
            )
            assert np.isclose(groups.contribution.sum(), row.itca)
        q = ns["audit"].table.query("group in ['Q0','Q1','Q2']")
        assert (q.contrast_contribution > q.reference_contribution).all()
        assert ns["audit"].metadata["screening_status"] == "disabled"
    finally:
        plt.close("all")
