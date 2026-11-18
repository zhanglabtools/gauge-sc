"""ITCA+ decisions preserve raw scores and explicit comparison scope."""

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from sklearn.tree import DecisionTreeClassifier

import gauge_sc as ip


@pytest.fixture
def evaluation():
    x = np.repeat(np.arange(3), 20)
    a = AnnData(
        x[:, None].astype(float),
        obs=pd.DataFrame(
            {"coarse": np.where(x == 0, "A", "B"), "fine": x.astype(str)},
            index=[f"c{i}" for i in range(60)],
        ),
    )
    r = ip.evaluate(
        a,
        groupby=["coarse", "fine"],
        views={"molecular": {"X": True}},
        classifier=DecisionTreeClassifier(random_state=0),
    )
    # Controlled scores isolate the decision rule from classifier behavior.
    r.table.loc[r.table["groupby"] == "coarse", "itca"] = 0.5
    r.table.loc[r.table["groupby"] == "fine", "itca"] = 0.58
    return r


def test_formula_path_no_refit_and_export(evaluation, monkeypatch, tmp_path):
    before = evaluation.table.copy(deep=True)

    def no_fit(*args, **kwargs):
        raise AssertionError("Selection must not retrain")

    monkeypatch.setattr(DecisionTreeClassifier, "fit", no_fit)
    selected = ip.select(evaluation, view="molecular", population_ratio=[0.01, 0.05])
    assert selected.path["groupby"].tolist() == ["fine", "coarse"]
    assert len(selected.table) == 4
    np.testing.assert_allclose(
        selected.table.delta,
        -selected.table.population_ratio * np.log(selected.table.population_ratio),
    )
    np.testing.assert_allclose(
        selected.table.itca_plus,
        selected.table.itca - selected.table.delta * selected.table.n_groups,
    )
    pd.testing.assert_frame_equal(evaluation.table, before)
    selected.save(tmp_path / "selection")
    saved = json.loads((tmp_path / "selection/result.json").read_text())
    assert saved["kind"] == "selection"
    assert saved["metadata"]["criterion"] == "ITCA+"
    assert len(saved["tables"]["path"]) == 2
    ax = ip.pl.table(selected, show=False)
    labels = [t.get_text() for t in ax.texts]
    assert {"ITCA", "ITCA+", "Penalty"} <= set(labels)
    plt.close("all")


def test_ties_and_monotone_complexity(evaluation):
    delta = -0.1 * np.log(0.1)
    evaluation.table["itca"] = delta * evaluation.table.n_groups
    tied = ip.select(evaluation, view="molecular", population_ratio=0.1)
    assert (tied.table.itca_plus == 0).all()
    assert tied.path.iloc[0]["groupby"] == "coarse"
    evaluation.table["itca"] = 0.5
    assert (
        ip.select(evaluation, view="molecular", population_ratio=0.01).path.iloc[0][
            "groupby"
        ]
        == "coarse"
    )
    evaluation.table["n_groups"] = 2
    a = ip.select(evaluation, view="molecular", population_ratio=0.01)
    evaluation.table = evaluation.table.iloc[::-1]
    b = ip.select(evaluation, view="molecular", population_ratio=0.01)
    assert a.path.iloc[0]["groupby"] == b.path.iloc[0]["groupby"] == "coarse"
    evaluation.table.loc[
        evaluation.table["groupby"] == "fine", ["n_groups", "itca"]
    ] = [3, 0.58]
    path = ip.select(
        evaluation, view="molecular", population_ratio=np.linspace(0.001, 1 / np.e, 20)
    ).path
    assert (np.diff(path.n_groups) <= 0).all()


@pytest.mark.parametrize(
    "p",
    [
        None,
        True,
        0,
        -0.1,
        0.5,
        float("nan"),
        float("inf"),
        [],
        [0.1, 0.1],
        [True, 0.1],
        [[0.1]],
        "0.1",
    ],
)
def test_invalid_population_ratio(evaluation, p):
    with pytest.raises(ValueError, match="POPULATION_RATIO"):
        ip.select(evaluation, view="molecular", population_ratio=p)


def test_selection_rejects_incomparable_or_invalid_candidates(evaluation):
    with pytest.raises(ValueError, match="MISSING_VIEW"):
        ip.select(evaluation, view="qc", population_ratio=0.1)
    mask = evaluation.predictions["groupby"] == "fine"
    evaluation.predictions.loc[mask, "obs_name"] += "_changed"
    with pytest.raises(ValueError, match="CELL_MISMATCH"):
        ip.select(evaluation, view="molecular", population_ratio=0.1)
    evaluation.table.loc[0, "itca"] = np.nan
    with pytest.raises(ValueError, match="INVALID_CANDIDATES"):
        ip.select(evaluation, view="molecular", population_ratio=0.1)


def test_direct_delta_matches_population_calibration(evaluation):
    p = 0.05
    calibrated = ip.select(evaluation, view="molecular", population_ratio=p)
    direct = ip.select(evaluation, view="molecular", delta=-p * np.log(p))
    np.testing.assert_allclose(calibrated.table.itca_plus, direct.table.itca_plus)
    assert calibrated.path["groupby"].tolist() == direct.path["groupby"].tolist()
    zero = ip.select(evaluation, view="molecular", delta=0)
    np.testing.assert_allclose(zero.table.itca_plus, zero.table.itca)
    ip.pl.table(direct, show=False)
    plt.close("all")
