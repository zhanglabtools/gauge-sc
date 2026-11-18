import json

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from sklearn.tree import DecisionTreeClassifier

import gauge_sc as ip
from gauge_sc._audit import _mirror_threshold

VIEWS = {"molecular": {"obsm": "X_pca"}, "qc": {"obs": ["depth"], "log1p": ["depth"]}}


@pytest.fixture
def data():
    labels = np.repeat(["A", "B", "C"], [31, 27, 22])
    codes = pd.Categorical(labels).codes.astype(float)
    obs = pd.DataFrame(
        {"labels": pd.Categorical(labels), "depth": np.ones(len(labels)) * 10},
        index=[f"cell-{i}" for i in range(len(labels))],
    )
    obj = ad.AnnData(sparse.csr_matrix(np.column_stack([codes, codes + 1])), obs=obs)
    obj.obsm["X_pca"] = codes[:, None]
    return obj


def run(data, **kwargs):
    return ip.evaluate(
        data,
        groupby="labels",
        views=VIEWS,
        classifier=DecisionTreeClassifier(random_state=0),
        **kwargs,
    )


def test_perfect_score_entropy_additivity_and_shared_folds(data):
    before = data.copy()
    result = run(data)
    molecular = result.table.query("view == 'molecular'").iloc[0]
    assert molecular.itca == pytest.approx(molecular.entropy)
    for row in result.table.itertuples():
        groups = result.groups.query("view == @row.view")
        assert groups.contribution.sum() == pytest.approx(row.itca)
        assert 0 <= row.itca <= row.entropy + 1e-12
    predictions = result.predictions.pivot(
        index="obs_name", columns="view", values="fold"
    )
    assert (predictions.molecular == predictions.qc).all()
    for view in VIEWS:
        p = result.predictions.query("view == @view")
        assert p.obs_name.is_unique and len(p) == data.n_obs
    pd.testing.assert_frame_equal(data.obs, before.obs)
    assert (data.X != before.X).nnz == 0
    assert data.uns == before.uns


def test_fold_mean_and_pooled_are_explicit(data):
    pooled = run(data, aggregation="pooled")
    mean = run(data, aggregation="fold_mean")
    np.testing.assert_allclose(pooled.groups.recall, pooled.groups.recall_pooled)
    np.testing.assert_allclose(mean.groups.recall, mean.groups.recall_fold_mean)
    assert pooled.metadata["aggregation"] == "pooled"
    qc = pooled.predictions.query("view == 'qc'")
    g = pooled.groups.query("view == 'qc'").set_index("group")
    for label in g.index:
        rows = qc[qc.group == label]
        assert g.loc[label, "recall"] == pytest.approx(
            (rows.group == rows.predicted_group).mean()
        )


def test_audit_formula_and_ineligible_groups(data):
    result = run(data)
    audit = ip.audit(result, reference="molecular", contrast="qc")
    assert len(audit.table) == 3
    c = audit.table.set_index("group").loc["C"]
    assert not c.eligible and pd.isna(c.contrast) and c.reason
    eligible = audit.table.query("eligible")
    difference = eligible.contrast_contribution - eligible.reference_contribution
    np.testing.assert_allclose(eligible.contrast, difference - difference.median())
    assert "candidate" not in audit.table
    assert audit.metadata["screening_status"] == "disabled"


def test_no_eligible_and_no_threshold(data):
    audit = ip.audit(
        run(data),
        reference="molecular",
        contrast="qc",
        min_test_cells=100,
        maximum_mirror_ratio=0.2,
    )
    assert audit.metadata["status"] == "no_eligible_groups"
    assert audit.metadata["threshold"] is None
    assert not audit.table.candidate.any()
    assert audit.table.contrast.isna().all()
    assert _mirror_threshold(np.array([0.0, 0.0]), 1) is None


def test_mirror_threshold_matches_manuscript_rule():
    values = np.array([-0.1, 0.0, 0.2, 0.3, 0.4])
    assert _mirror_threshold(values, 1 / 3) == pytest.approx(0.2)
    assert _mirror_threshold(values, 0.2) is None


def test_writeback_reorder_h5ad_and_collision(data, tmp_path):
    audit = ip.audit(
        run(data), reference="molecular", contrast="qc", maximum_mirror_ratio=1
    )
    shuffled = data[::-1].copy()
    audit.write_to_adata(shuffled)
    lookup = audit.table.set_index("group").contrast
    expected = shuffled.obs.labels.astype(str).map(lookup)
    np.testing.assert_allclose(
        shuffled.obs.gauge_sc_audit_contrast, expected, equal_nan=True
    )
    path = tmp_path / "audit.h5ad"
    shuffled.write_h5ad(path)
    restored = ad.read_h5ad(path)
    assert restored.uns["gauge_sc_audit"]["score_level"].startswith("group")
    with pytest.raises(ValueError, match="OUTPUT_EXISTS"):
        audit.write_to_adata(shuffled)
    changed = data.copy()
    changed.obs["labels"] = "different"
    with pytest.raises(ValueError, match="GROUPING_MISMATCH"):
        audit.write_to_adata(changed)
    with pytest.raises(ValueError, match="CELL_MISMATCH"):
        audit.write_to_adata(data[:-1].copy())


def test_export_and_plots(data, tmp_path):
    result = run(data)
    audit = ip.audit(result, reference="molecular", contrast="qc")
    for obj, name in [(result, "evaluation"), (audit, "audit")]:
        obj.save(tmp_path / name)
        loaded = json.loads((tmp_path / name / "result.json").read_text())
        assert loaded["kind"] == name
        assert loaded["schema_version"] == 1
        assert len(loaded["tables"]["table"]) == len(obj.table)
        with pytest.raises(FileExistsError):
            obj.save(tmp_path / name)
        ax = ip.pl.table(obj, show=False)
        ax.figure.savefig(tmp_path / f"{name}-table.svg")
    fig, axes = plt.subplots(1, 2)
    assert ip.pl.audit(audit, ax=axes, show=False) is fig
    fig.savefig(tmp_path / "audit.svg")
    assert (tmp_path / "audit.svg").stat().st_size > 1000
    plt.close("all")


def test_multiple_candidates_require_explicit_audit(data):
    data.obs["other"] = data.obs.labels
    result = ip.evaluate(
        data,
        groupby=["labels", "other"],
        views=VIEWS,
        classifier=DecisionTreeClassifier(random_state=0),
    )
    assert len(result.table) == 4
    with pytest.raises(ValueError, match="AMBIGUOUS_GROUPING"):
        ip.audit(result, reference="molecular", contrast="qc")
    assert (
        len(
            ip.audit(
                result, reference="molecular", contrast="qc", groupby="other"
            ).table
        )
        == 3
    )


@pytest.mark.parametrize(
    "case, message",
    [
        ("tiny", "INSUFFICIENT_GROUP_SIZE"),
        ("missing", "MISSING_LABEL"),
        ("one", "SINGLE_GROUP"),
        ("nonfinite", "NONFINITE_VIEW"),
        ("negative", "INVALID_TRANSFORM"),
    ],
)
def test_invalid_inputs_are_explained(data, case, message):
    obs = data.obs
    if case == "tiny":
        obs["labels"] = ["tiny"] + list(obs.labels)[1:]
    elif case == "missing":
        obs["labels"] = [None] + list(obs.labels)[1:]
    elif case == "one":
        obs["labels"] = "one"
    elif case == "nonfinite":
        obs["depth"] = np.nan
    else:
        obs["depth"] = -1
    with pytest.raises(ValueError, match=message):
        run(data)


def test_sparse_X_and_layer_views(data):
    data.layers["copy"] = data.X.copy()
    result = ip.evaluate(
        data,
        groupby="labels",
        views={"x": {"X": True}, "layer": {"layer": "copy"}},
        classifier=DecisionTreeClassifier(random_state=0),
    )
    assert result.table.itca.iloc[0] == result.table.itca.iloc[1]


def test_example_is_deterministic_and_synthetic():
    a, b = ip.datasets.audit_example(), ip.datasets.audit_example()
    assert a.shape == (560, 50)
    assert a.uns["gauge_sc_example"]["kind"] == "synthetic"
    assert (a.X != b.X).nnz == 0
    pd.testing.assert_frame_equal(a.obs, b.obs)
    np.testing.assert_array_equal(a.obsm["X_pca"], b.obsm["X_pca"])
    assert "counts" not in a.layers
    assert "not derived from X" in a.uns["gauge_sc_example"]["qc_provenance"]


def test_reproducibility_and_feature_provenance(data):
    first, second = run(data), run(data)
    pd.testing.assert_frame_equal(first.predictions, second.predictions)
    assert first.metadata["view_fingerprints"] == second.metadata["view_fingerprints"]
    data.obsm["X_pca"] = data.obsm["X_pca"] + 0.1
    third = run(data)
    assert (
        first.metadata["view_fingerprints"]["molecular"]
        != third.metadata["view_fingerprints"]["molecular"]
    )


def test_random_forest_default_is_runnable():
    data = ip.datasets.audit_example()
    result = ip.evaluate(data, groupby="group", views={"molecular": {"obsm": "X_pca"}})
    assert result.table.itca.iloc[0] > 0


@pytest.mark.parametrize("seed", range(5))
def test_teaching_example_screens_designed_qc_groups(seed):
    data = ip.datasets.audit_example(random_state=seed)
    result = ip.evaluate(
        data,
        groupby="group",
        views={
            "molecular": {"obsm": "X_pca"},
            "qc": {
                "obs": ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
                "log1p": ["total_counts", "n_genes_by_counts"],
            },
        },
    )
    audit = ip.audit(
        result, reference="molecular", contrast="qc", maximum_mirror_ratio=1 / 3
    )
    assert audit.table.eligible.all()
    assert audit.metadata["screening_status"] == "threshold_found"
    expected = set(data.obs.loc[data.obs.simulation_role == "qc_associated", "group"])
    assert expected <= set(audit.table.loc[audit.table.candidate, "group"])
    assert audit.table.loc[audit.table.candidate, "contrast"].gt(0).all()
    # Screening is still opt-in; the tutorial does not alter the API default.
    unscreened = ip.audit(result, reference="molecular", contrast="qc")
    assert "candidate" not in unscreened.table


@pytest.mark.parametrize(
    "values, ratio, expected_reason",
    [
        ([-0.1, 0, 0.2, 0.3, 0.4], None, "disabled"),
        ([-0.1, 0, 0.2, 0.3, 0.4], 1 / 3, "threshold_found"),
        ([-0.1, 0, 0.2, 0.3, 0.4], 0.2, "insufficient_positive_groups"),
        ([-0.4, -0.3, -0.2, 0, 0.2, 0.3, 0.4], 1 / 3, "mirror_ratio_not_met"),
        ([], 1 / 3, "no_eligible_groups"),
        ([-0.1, 0, 0.2, 0.3, 0.4], 0.33, "insufficient_positive_groups"),
    ],
)
def test_screening_diagnostics(values, ratio, expected_reason):
    from gauge_sc._audit import _screening_summary

    values = np.asarray(values)
    threshold = None if ratio is None else _mirror_threshold(values, ratio)
    detail = _screening_summary(values, ratio, threshold)
    assert detail["reason"] == expected_reason
    assert detail["positive_group_count"] == (values > 0).sum()
    json.dumps(detail, allow_nan=False)
    if threshold is None:
        assert detail["positive_tail_count"] is None
        assert detail["negative_tail_count"] is None
        assert detail["achieved_mirror_ratio"] is None
    else:
        assert detail["positive_tail_count"] == (values >= threshold).sum()
        assert detail["negative_tail_count"] == (values <= -threshold).sum()
        assert detail["achieved_mirror_ratio"] <= ratio
    if ratio == 1 / 3:
        assert detail["minimum_positive_count"] == 3
    if ratio == 0.33:
        assert detail["minimum_positive_count"] == 4


def test_sensitivity_reuses_scores_and_preserves_serialized_explanation(
    data, tmp_path, monkeypatch
):
    result = run(data)
    before = result.groups.copy(deep=True)

    def cannot_fit(*args, **kwargs):
        raise AssertionError("Sensitivity analysis must not refit a classifier")

    monkeypatch.setattr(DecisionTreeClassifier, "fit", cannot_fit)
    audits = [
        ip.audit(result, reference="molecular", contrast="qc", maximum_mirror_ratio=q)
        for q in [0.2, 1 / 3, 0.5]
    ]
    for audit in audits:
        np.testing.assert_allclose(
            audit.table.contrast, audits[0].table.contrast, equal_nan=True
        )
    pd.testing.assert_frame_equal(result.groups, before)
    audit = audits[1]
    audit.save(tmp_path / "diagnostics")
    saved = json.loads((tmp_path / "diagnostics/result.json").read_text())
    assert saved["metadata"]["screening"] == audit.metadata["screening"]
    audit.write_to_adata(data)
    data.write_h5ad(tmp_path / "diagnostics.h5ad")
    restored = ad.read_h5ad(tmp_path / "diagnostics.h5ad")
    assert (
        json.loads(restored.uns["gauge_sc_audit"]["metadata_json"])["screening"]
        == audit.metadata["screening"]
    )
