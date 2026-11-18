"""Plot controls preserve scientific output and embedding behavior."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import gauge_sc as ip
from gauge_sc._results import AuditResult


@pytest.fixture
def result():
    return AuditResult(
        table=pd.DataFrame(
            {
                "group": ["A", "B", "C", "D"],
                "n_cells": [100, 80, 60, 10],
                "reference_contribution": [0.2, 0.21, 0.1, 0.05],
                "contrast_contribution": [0.21, 0.22, 0.3, 0.1],
                "contrast": [-0.1, 0.0, 0.15, np.nan],
                "eligible": [True, True, True, False],
                "candidate": [False, False, True, False],
            }
        ),
        assignments=pd.DataFrame(),
        metadata={"reference": "molecular", "contrast": "qc", "threshold": 0.12},
    )


@pytest.fixture(autouse=True)
def cleanup():
    yield
    plt.close("all")


def test_panels_styles_and_no_mutation(result):
    original = result.table.copy(deep=True)
    rc = matplotlib.rcParams.copy()
    for panel, count in [("both", 2), ("contributions", 1), ("contrast", 1)]:
        fig = ip.pl.audit(
            result,
            panel=panel,
            figsize=(8, 4),
            show=False,
            scatter_kws={"s": 72, "alpha": 0.6},
        )
        assert len(fig.axes) == count
        assert all(ax.figure is fig for ax in fig.axes)
        for ax in fig.axes:
            for points in ax.collections:
                assert list(points.get_sizes()) == [72]
                assert points.get_alpha() == 0.6
    pd.testing.assert_frame_equal(original, result.table)
    for key in rc:
        if key != "backend":
            assert matplotlib.rcParams[key] == rc[key], key


def test_candidates_explicit_labels_and_input_order(result):
    fig = ip.pl.audit(result, show=False)
    assert [t.get_text() for t in fig.axes[0].texts] == ["C"]
    assert [t.get_text() for t in fig.axes[1].get_yticklabels()] == ["C", "B", "A", "D"]
    fig = ip.pl.audit(
        result, annotate=["A", "B"], sort="input", text_kws={"fontsize": 13}, show=False
    )
    assert [t.get_text() for t in fig.axes[0].texts] == ["A", "B"]
    assert all(t.get_fontsize() == 13 for t in fig.axes[0].texts)
    assert [t.get_text() for t in fig.axes[1].get_yticklabels()] == ["A", "B", "C", "D"]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    a, b = fig.axes[0].texts
    assert not a.get_window_extent(renderer).overlaps(b.get_window_extent(renderer))


def test_single_panel_axes_and_conflicts(result):
    fig, ax = plt.subplots()
    assert ip.pl.audit(result, panel="contrast", ax=ax, show=False) is fig
    with pytest.raises(ValueError, match="figsize"):
        ip.pl.audit(result, panel="contrast", ax=ax, figsize=(4, 4))
    with pytest.raises(ValueError, match="requires 2"):
        ip.pl.audit(result, ax=ax)
    with pytest.raises(ValueError, match="distinct axes"):
        ip.pl.audit(result, ax=[ax, ax])
    with pytest.raises(ValueError, match="figsize"):
        ip.pl.table(result, ax=ax, figsize=(4, 4))


@pytest.mark.parametrize(
    "kwargs",
    [
        {"panel": "bad"},
        {"annotate": "bad"},
        {"sort": "bad"},
        {"annotate": ["unknown"]},
        {"scatter_kws": {"colour": "red"}},
        {"scatter_kws": {"c": "red"}},
        {"text_kws": {"xytext": (0, 0)}},
        {"scatter_kws": {"s": [10, 20]}},
        {"scatter_kws": {"s": "large"}},
    ],
)
def test_invalid_options_fail_before_plotting(result, kwargs):
    before = plt.get_fignums()
    with pytest.raises(ValueError):
        ip.pl.audit(result, **kwargs)
    assert plt.get_fignums() == before


def test_missing_values_not_shown_as_zero_and_table_status(result):
    fig = ip.pl.audit(result, panel="contrast", show=False)
    ax = fig.axes[0]
    assert sum(len(c.get_offsets()) for c in ax.collections) == 3
    assert "Not eligible" in [t.get_text() for t in ax.texts]
    table = ip.pl.table(result, show=False)
    text = [t.get_text() for t in table.texts]
    assert "—" in text and "Candidate" in text and "Not eligible" in text
    result.table = result.table.drop(columns="candidate")
    table = ip.pl.table(result, show=False)
    assert "Not screened" in [t.get_text() for t in table.texts]


def test_palette_mapping_stays_aligned_after_sorting(result):
    from matplotlib.colors import to_rgba

    palette = {"A": "#1f77b4", "B": "#ff7f0e", "C": "#279e68", "D": "#d62728"}
    fig = ip.pl.audit(result, palette=palette, show=False)
    np.testing.assert_allclose(
        fig.axes[0].collections[0].get_facecolors()[0][:3], to_rgba(palette["A"])[:3]
    )
    ordered = ["C", "B", "A"]
    for points, group in zip(fig.axes[1].collections, ordered):
        np.testing.assert_allclose(
            points.get_facecolors()[0][:3], to_rgba(palette[group])[:3]
        )
    assert fig.axes[1].get_yticklabels()[0].get_fontweight() == "bold"
    ax = ip.pl.table(result, palette=palette, show=False)
    labels = {t.get_text(): t.get_color() for t in ax.texts if t.get_text() in palette}
    assert {to_rgba(value) for value in labels.values()} == {to_rgba("black")}
    swatches = [p for p in ax.patches if p.get_width() < 0.01]
    for patch, expected in zip(swatches, palette.values()):
        np.testing.assert_allclose(patch.get_facecolor(), to_rgba(expected))


def test_palette_validation_and_monochrome_override(result):
    with pytest.raises(ValueError, match="missing groups"):
        ip.pl.audit(result, palette={"A": "red"}, show=False)
    with pytest.raises(ValueError, match="at least one color"):
        ip.pl.audit(result, palette=["red"], show=False)
    fig = ip.pl.audit(result, color="black", show=False)
    np.testing.assert_allclose(fig.axes[0].collections[0].get_facecolors()[:3, :3], 0)


def test_named_context_and_font_scaling(result):
    import seaborn as sns

    expected = sns.plotting_context("paper", font_scale=1.25)
    fig = ip.pl.audit(
        result,
        context="paper",
        font_scale=1.25,
        fontfamily="DejaVu Sans Mono",
        annotate="all",
        show=False,
    )
    ax = fig.axes[0]
    assert ax.xaxis.label.get_fontsize() == pytest.approx(expected["axes.labelsize"])
    assert ax.spines["left"].get_linewidth() == pytest.approx(
        expected["axes.linewidth"]
    )
    assert ax.get_xticklabels()[0].get_fontsize() == pytest.approx(
        expected["xtick.labelsize"]
    )
    assert ax.get_xticklabels()[0].get_fontfamily() == ["DejaVu Sans Mono"]
    assert ax.texts[0].get_fontfamily() == ["DejaVu Sans Mono"]
    ax = ip.pl.table(
        result,
        context="paper",
        font_scale=1.25,
        fontfamily="DejaVu Sans Mono",
        show=False,
    )
    body = next(t for t in ax.texts if t.get_text() == "A")
    assert body.get_fontsize() == pytest.approx(expected["font.size"])
    assert body.get_fontfamily() == ["DejaVu Sans Mono"]


def test_inherit_external_context_and_restore_on_error(result):
    import seaborn as sns

    with (
        sns.plotting_context("talk"),
        matplotlib.rc_context(
            {
                "axes.edgecolor": "purple",
                "axes.labelcolor": "navy",
                "text.color": "navy",
                "axes.grid": True,
                "axes.spines.top": True,
                "font.family": ["DejaVu Serif"],
            }
        ),
    ):
        before = matplotlib.rcParams.copy()
        original, provided = plt.subplots()
        fig = ip.pl.audit(
            result, panel="contrast", ax=provided, context=None, show=False
        )
        assert fig is original
        assert provided.xaxis.label.get_fontsize() == before["axes.labelsize"]
        assert provided.xaxis.label.get_fontfamily() == ["DejaVu Serif"]
        assert provided.spines["top"].get_visible()
        assert any(g.get_visible() for g in provided.get_xgridlines())
        assert provided.spines["left"].get_edgecolor() == matplotlib.colors.to_rgba(
            "purple"
        )
        with pytest.raises(ValueError):
            ip.pl.audit(result, context="paper", font_scale=2, panel="invalid")
        for key in before:
            assert matplotlib.rcParams[key] == before[key], key


@pytest.mark.parametrize(
    "kwargs",
    [
        {"font_scale": 0},
        {"font_scale": float("nan")},
        {"context": "invalid"},
        {"fontfamily": ""},
    ],
)
def test_invalid_context_parameters(result, kwargs):
    with pytest.raises(ValueError):
        ip.pl.audit(result, show=False, **kwargs)


def test_many_groups_facet_without_dropping_or_recomputing(result):
    result.table = pd.concat([result.table] * 6, ignore_index=True)
    result.table["group"] = [f"G{i}" for i in range(24)]
    before = result.table.copy(deep=True)
    fig = ip.pl.audit(result, palette=None, show=False)
    assert len(fig.axes) == 4
    contrasts = fig.axes[1:]
    expected = result.table.sort_values(
        "contrast", ascending=False, na_position="last", kind="stable"
    ).group.tolist()
    assert [t.get_text() for ax in contrasts for t in ax.get_yticklabels()] == expected
    assert all(ax.get_xlim() == contrasts[0].get_xlim() for ax in contrasts)
    assert all(ax.get_ylim() == contrasts[0].get_ylim() for ax in contrasts)
    fig.canvas.draw()
    heights = [ax.get_window_extent().height for ax in fig.axes]
    assert max(heights) - min(heights) < 1
    pd.testing.assert_frame_equal(result.table, before)


@pytest.mark.parametrize("value", [0, -1, True, 2.5])
def test_invalid_groups_per_panel(result, value):
    with pytest.raises(ValueError, match="groups_per_panel"):
        ip.pl.audit(result, groups_per_panel=value, show=False)
