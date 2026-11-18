"""Locally styled Matplotlib views of stored ITCA+ results.

Plotting never refits models or changes screening. Local rcParams changes are
restored on return, including when plotting raises an exception.
"""

from collections.abc import Mapping
from functools import wraps
from numbers import Integral, Real
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, fontManager
import seaborn as sns
from matplotlib.axes import Axes
from matplotlib.colors import to_rgba
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

_INK = "#000000"
_COLOR = "#487D8B"
_CANDIDATE = "#B85C70"


# Scanpy's reordered/adjusted Vega20 categorical palette. Kept as color values
# so core plotting does not require installing Scanpy.
_VEGA20 = [
    "#1f77b4",
    "#ff7f0e",
    "#279e68",
    "#d62728",
    "#aa40fc",
    "#8c564b",
    "#e377c2",
    "#b5bd61",
    "#17becf",
    "#aec7e8",
    "#ffbb78",
    "#98df8a",
    "#ff9896",
    "#c5b0d5",
    "#c49c94",
    "#f7b6d2",
    "#dbdb8d",
    "#9edae5",
    "#ad494a",
    "#8c6d31",
]


_INTER_REGISTERED = False


def _register_inter():
    """Register bundled Inter faces only when explicitly requested."""
    global _INTER_REGISTERED
    if not _INTER_REGISTERED:
        for path in sorted((Path(__file__).parent / "fonts").glob("*.otf")):
            fontManager.addfont(str(path))
        _INTER_REGISTERED = True


def _points(value):
    return FontProperties(size=value).get_size_in_points()


def _with_context(function):
    """Apply an exception-safe local context; keep the public signature intact."""

    @wraps(function)
    def wrapped(*args, **kwargs):
        context = kwargs.get("context", "paper")
        scale = kwargs.get("font_scale", 1)
        family = kwargs.get("fontfamily")
        if (
            isinstance(scale, bool)
            or not isinstance(scale, Real)
            or not np.isfinite(scale)
            or scale <= 0
        ):
            raise ValueError("font_scale must be a positive finite number.")
        if context not in {None, "paper", "notebook", "talk", "poster"}:
            raise ValueError("context must be paper, notebook, talk, poster, or None.")
        if family is not None and (not isinstance(family, str) or not family.strip()):
            raise ValueError("fontfamily must be a nonempty font family name or None.")
        if context is None:
            rc = {}
            for key in [
                "font.size",
                "axes.titlesize",
                "axes.labelsize",
                "xtick.labelsize",
                "ytick.labelsize",
                "legend.fontsize",
                "legend.title_fontsize",
            ]:
                rc[key] = _points(mpl.rcParams[key]) * scale
        else:
            rc = dict(sns.plotting_context(context, font_scale=scale))
            rc.update(
                {
                    "axes.edgecolor": "black",
                    "axes.labelcolor": "black",
                    "text.color": "black",
                    "xtick.color": "black",
                    "ytick.color": "black",
                    "xtick.labelcolor": "black",
                    "ytick.labelcolor": "black",
                    "axes.facecolor": "white",
                    "figure.facecolor": "white",
                    "axes.grid": False,
                    "axes.spines.top": False,
                    "axes.spines.right": False,
                    "axes.spines.bottom": True,
                    "axes.spines.left": True,
                    "xtick.direction": "out",
                    "ytick.direction": "out",
                }
            )
        if family == "Inter":
            _register_inter()
        if family is not None:
            rc["font.family"] = [family]
        with mpl.rc_context(rc=rc):
            return function(*args, **kwargs)

    return wrapped


def _group_colors(groups, palette, color=None):
    """Resolve colors before sorting; explicit mappings preserve external order."""
    if color is not None or palette is None:
        chosen = color if color is not None else _COLOR
        to_rgba(chosen)
        return dict.fromkeys(groups, chosen)
    if isinstance(palette, str):
        if palette != "vega20":
            raise ValueError(
                "palette must be vega20, a color sequence, a group-color mapping, or None."
            )
        palette = _VEGA20
    if isinstance(palette, Mapping):
        missing = set(groups) - set(palette)
        if missing:
            raise ValueError(f"palette is missing groups: {sorted(missing)}")
        resolved = {group: palette[group] for group in groups}
    else:
        if not isinstance(palette, (list, tuple, np.ndarray)) or len(palette) < len(
            groups
        ):
            raise ValueError(
                "palette needs at least one color per group; supply a longer palette or explicit mapping."
            )
        resolved = dict(zip(groups, palette))
    for value in resolved.values():
        to_rgba(value)
    return resolved


def _view_label(name):
    return "QC" if name.lower() == "qc" else name.replace("_", " ").capitalize()


def _style_options(options, allowed, name):
    if options is None:
        return {}
    if not isinstance(options, Mapping):
        raise ValueError(f"{name} must be a dictionary of styling parameters.")
    unknown = set(options) - set(allowed)
    if unknown:
        raise ValueError(
            f"Unsupported {name} keys: {sorted(map(str, unknown))}. "
            f"Allowed: {sorted(allowed)}. Use explicit color/layout parameters "
            "or edit the returned Matplotlib artists."
        )
    return dict(options)


def _axes(ax, panel, figsize, rows, groups_per_panel):
    contrasts = (rows + groups_per_panel - 1) // groups_per_panel
    count = (1 if panel != "contrast" else 0) + (
        contrasts if panel != "contributions" else 0
    )
    if ax is not None:
        if figsize is not None:
            raise ValueError(
                "figsize cannot be used with ax; size the supplied figure."
            )
        axes = np.asarray([ax] if isinstance(ax, Axes) else ax, dtype=object).ravel()
        if len(axes) != count or not all(isinstance(a, Axes) for a in axes):
            raise ValueError(f"panel={panel!r} requires {count} Matplotlib axes.")
        if len({id(a) for a in axes}) != count or any(
            a.figure is not axes[0].figure for a in axes
        ):
            raise ValueError("Provide distinct axes on the same figure.")
        return axes[0].figure, axes
    columns = 2 if count == 4 else min(count, 3)
    rows = (count + columns - 1) // columns
    fig = plt.figure(
        figsize=figsize or (3.6 * columns, 3.6 * rows), layout="constrained"
    )
    grid = fig.add_gridspec(rows, columns)
    axes = np.asarray(
        [fig.add_subplot(grid[i // columns, i % columns]) for i in range(count)]
    )
    fig.set_facecolor(mpl.rcParams["figure.facecolor"])
    return fig, axes


def _style_axis(ax):
    rc = mpl.rcParams
    ax.set_facecolor(rc["axes.facecolor"])
    ax.set_axisbelow(True)
    for side in ["left", "right", "top", "bottom"]:
        ax.spines[side].set_visible(rc[f"axes.spines.{side}"])
        ax.spines[side].set_color(rc["axes.edgecolor"])
        ax.spines[side].set_linewidth(rc["axes.linewidth"])
    for which in ["x", "y"]:
        prefix = "xtick" if which == "x" else "ytick"
        ax.tick_params(
            axis=which,
            colors=rc[f"{prefix}.color"],
            labelcolor=(
                rc[f"{prefix}.color"]
                if rc[f"{prefix}.labelcolor"] == "inherit"
                else rc[f"{prefix}.labelcolor"]
            ),
            labelsize=rc[f"{prefix}.labelsize"],
            length=rc[f"{prefix}.major.size"],
            width=rc[f"{prefix}.major.width"],
            direction=rc[f"{prefix}.direction"],
        )
        axis = ax.xaxis if which == "x" else ax.yaxis
        axis.label.set(
            color=rc["axes.labelcolor"],
            fontsize=rc["axes.labelsize"],
            fontfamily=rc["font.family"],
        )
    ax.grid(False, which="both")
    if rc["axes.grid"]:
        ax.grid(
            True,
            axis=rc["axes.grid.axis"],
            which=rc["axes.grid.which"],
            color=rc["grid.color"],
            linewidth=rc["grid.linewidth"],
        )


def _tick_fonts(ax):
    # Also update supplied Axes, whose tick Text objects predate this context.
    for label in [
        *ax.get_xticklabels(),
        *ax.get_yticklabels(),
        ax.xaxis.get_offset_text(),
        ax.yaxis.get_offset_text(),
    ]:
        label.set_fontfamily(mpl.rcParams["font.family"])


def _labels(ax, frame, labels, text_kws):
    """Place selected labels, shifting colliding text boxes in display space."""
    placed = []
    for row in frame.itertuples():
        if row.group not in labels:
            continue
        annotation = ax.annotate(
            row.group,
            (row.reference_contribution, row.contrast_contribution),
            xytext=(7, 7),
            textcoords="offset points",
            color=mpl.rcParams["text.color"],
            annotation_clip=False,
            **text_kws,
        )
        placed.append(annotation)
    if placed:
        ax.figure.canvas.draw()
        renderer = ax.figure.canvas.get_renderer()
        for i, annotation in enumerate(placed):
            for dx in (7, 30, 53, 76, -30, -53):
                found = False
                for dy in (7, 21, 35, -14, -28):
                    annotation.set_position((dx, dy))
                    box = annotation.get_window_extent(renderer).expanded(1.12, 1.2)
                    bounds = ax.get_window_extent(renderer)
                    if (
                        box.x0 >= bounds.x0
                        and box.x1 <= bounds.x1
                        and box.y0 >= bounds.y0
                        and box.y1 <= bounds.y1
                        and not any(
                            box.overlaps(other.get_window_extent(renderer))
                            for other in placed[:i]
                        )
                    ):
                        found = True
                        break
                if found:
                    break


@_with_context
def audit(
    result,
    *,
    panel="both",
    groups_per_panel=10,
    annotate="candidates",
    sort="contrast",
    figsize=None,
    palette="vega20",
    color=None,
    candidate_color=_CANDIDATE,
    scatter_kws=None,
    text_kws=None,
    ax=None,
    show=True,
    context="paper",
    font_scale=1,
    fontfamily=None,
):
    """Plot stored view contributions and centered contrasts; return Figure.

    Args:
        result (AuditResult): Existing audit. Screening is never recomputed.
        panel (str): 'both', 'contributions', or 'contrast'.
        groups_per_panel (int): Maximum rows per contrast facet (default 10).
            All groups remain visible; facets share limits and wrap at three columns.
        annotate (str or list[str]): Scatter labels: 'candidates' (default),
            'all', 'none', or explicit group names. Unscreened results have no
            candidate labels. Contrast rows always show group names.
        sort (str): Contrast rows in descending 'contrast' order, or input order.
        figsize (tuple, optional): Figure size in inches. Cannot accompany ax.
        palette (str, sequence or mapping): Scanpy-adjusted 'vega20' by default.
            Colors identify groups, resolved in input order before sorting.
            Use a group-color mapping to match an existing Scanpy palette.
            None uses a single color; insufficient palettes raise an error.
        color (color, optional): Explicit monochrome override of palette.
        candidate_color (color): Screening-threshold line color. Candidates also
            have bold contrast-row labels; their group colors remain unchanged.
        scatter_kws (dict, optional): Point styling: s (one constant area), alpha,
            marker, linewidths, edgecolors. Applies to both panels.
        text_kws (dict, optional): Scatter-label styling: fontsize, fontweight,
            fontstyle, fontfamily, alpha. Coordinates and colors are reserved.
        ax (Axes or sequence[Axes], optional): One axis for a single panel or
            one per contrast facet plus one for contributions when requested.
        show (bool): Whether to display the figure. Defaults to True.
        context (str or None): Seaborn paper (default), notebook, talk or poster.
            None inherits the active Matplotlib/Seaborn context and axis style.
        font_scale (float): Positive relative font multiplier; lines are unchanged.
        fontfamily (str, optional): Unified font family; None inherits the environment.
            Per-annotation text_kws values and supplied fontsize values override defaults.

    Notes:
        Positive centered contrasts are relative to the eligible-group median;
        they need not imply positive raw differences. Ineligible rows are shown
        as unavailable, never plotted as zero. Colors do not encode good/bad.
    """
    if panel not in {"both", "contributions", "contrast"}:
        raise ValueError("panel must be both, contributions, or contrast.")
    if (
        isinstance(groups_per_panel, bool)
        or not isinstance(groups_per_panel, (int, np.integer))
        or groups_per_panel < 1
    ):
        raise ValueError("groups_per_panel must be a positive integer.")
    if sort not in {"contrast", "input"}:
        raise ValueError("sort must be contrast or input.")
    scatter = {
        "s": mpl.rcParams["lines.markersize"] ** 2,
        "alpha": 1.0,
        "marker": "o",
        "linewidths": mpl.rcParams["patch.linewidth"],
        "edgecolors": "white",
    }
    scatter.update(
        _style_options(
            scatter_kws,
            {"s", "alpha", "marker", "linewidths", "edgecolors"},
            "scatter_kws",
        )
    )
    if (
        isinstance(scatter["s"], bool)
        or not isinstance(scatter["s"], Real)
        or not np.isfinite(scatter["s"])
        or scatter["s"] <= 0
    ):
        raise ValueError("scatter_kws['s'] must be one positive finite marker area.")
    text = {
        "fontsize": mpl.rcParams["xtick.labelsize"],
        "fontfamily": mpl.rcParams["font.family"],
    }
    text.update(
        _style_options(
            text_kws,
            {"fontsize", "fontweight", "fontstyle", "fontfamily", "alpha"},
            "text_kws",
        )
    )
    to_rgba(candidate_color)
    frame = result.table.copy()
    group_colors = _group_colors(list(frame.group), palette, color)
    if isinstance(annotate, str):
        if annotate not in {"candidates", "all", "none"}:
            raise ValueError(
                "annotate must be candidates, all, none, or a list of group names."
            )
        labels = (
            set(frame.group)
            if annotate == "all"
            else (
                set(
                    frame.loc[
                        frame.get("candidate", pd.Series(False, index=frame.index)),
                        "group",
                    ]
                )
                if annotate == "candidates"
                else set()
            )
        )
    elif isinstance(annotate, (list, tuple)) and all(
        isinstance(x, str) for x in annotate
    ):
        labels = set(annotate)
        if labels - set(frame.group):
            raise ValueError(
                f"Unknown annotation groups: {sorted(labels - set(frame.group))}"
            )
    else:
        raise ValueError("annotate must be a mode or list of group names.")
    fig, axes = _axes(ax, panel, figsize, len(frame), groups_per_panel)
    reference = _view_label(result.metadata["reference"])
    compared = _view_label(result.metadata["contrast"])
    for axis in axes:
        _style_axis(axis)
    if panel in {"both", "contributions"}:
        left = axes[0]
        colors = [
            group_colors[row.group] if row.eligible else "#B8C1C5"
            for row in frame.itertuples()
        ]
        limit = (
            max(
                frame.reference_contribution.max(),
                frame.contrast_contribution.max(),
                0.01,
            )
            * 1.18
        )
        left.plot(
            [0, limit],
            [0, limit],
            color=mpl.rcParams["axes.edgecolor"],
            lw=mpl.rcParams["axes.linewidth"] * 0.6,
            linestyle=(0, (3, 4)),
            zorder=1,
        )
        left.scatter(
            frame.reference_contribution,
            frame.contrast_contribution,
            c=colors,
            zorder=3,
            **scatter,
        )
        left.set(
            xlim=(0, limit),
            ylim=(0, limit),
            xlabel=f"{reference} contribution (nats)",
            ylabel=f"{compared} contribution (nats)",
        )
        left.set_aspect("equal", adjustable="box")
        left.xaxis.set_major_locator(MaxNLocator(4))
        left.yaxis.set_major_locator(MaxNLocator(4))
        left.set_title(
            "View contributions",
            loc="left",
            fontsize=mpl.rcParams["axes.titlesize"],
            fontfamily=mpl.rcParams["font.family"],
            color=mpl.rcParams["text.color"],
            fontweight=mpl.rcParams["axes.titleweight"],
            pad=10,
        )
    if panel in {"both", "contrast"}:
        ordered = (
            frame.sort_values(
                "contrast", ascending=False, na_position="last", kind="stable"
            )
            if sort == "contrast"
            else frame
        )
        ordered = ordered.reset_index(drop=True)
        finite = ordered.loc[ordered.eligible, "contrast"].to_numpy()
        threshold = result.metadata["threshold"]
        extent = (
            max(
                float(np.max(np.abs(finite))) if len(finite) else 0,
                threshold or 0,
                0.001,
            )
            * 1.22
        )
        contrast_axes = axes[1:] if panel == "both" else axes
        for facet, right in enumerate(contrast_axes):
            start = facet * groups_per_panel
            chunk = ordered.iloc[start : start + groups_per_panel].reset_index(
                drop=True
            )
            for i, row in chunk.iterrows():
                if row.eligible:
                    c = group_colors[row.group]
                    right.plot(
                        [0, row.contrast],
                        [i, i],
                        color=c,
                        lw=mpl.rcParams["lines.linewidth"],
                        zorder=2,
                    )
                    right.scatter([row.contrast], [i], c=[c], zorder=3, **scatter)
                else:
                    right.text(
                        0.98,
                        i,
                        "Not eligible",
                        transform=right.get_yaxis_transform(),
                        ha="right",
                        va="center",
                        fontsize=mpl.rcParams["ytick.labelsize"],
                        color=mpl.rcParams["text.color"],
                    )
            right.axvline(
                0,
                color=mpl.rcParams["axes.edgecolor"],
                lw=mpl.rcParams["axes.linewidth"] * 0.6,
                linestyle=(0, (3, 3)),
                zorder=1,
            )
            if threshold is not None:
                right.axvline(
                    threshold,
                    color=candidate_color,
                    lw=mpl.rcParams["lines.linewidth"],
                    linestyle=(0, (3, 3)),
                    label="Screening threshold",
                )
            right.set(
                xlim=(-extent, extent),
                ylim=(
                    (
                        groups_per_panel - 0.5
                        if len(ordered) > groups_per_panel
                        else len(chunk) - 0.5
                    ),
                    -0.5,
                ),
                yticks=np.arange(len(chunk)),
                yticklabels=chunk.group,
                xlabel=f"{compared} − {reference.lower()} (nats)\nMedian-centered",
            )
            right.set_title(
                (
                    "Group contrast"
                    if len(ordered) <= groups_per_panel
                    else f"Group contrast ({start + 1}–{start + len(chunk)})"
                ),
                loc="left",
                fontsize=mpl.rcParams["axes.titlesize"],
                fontfamily=mpl.rcParams["font.family"],
                color=mpl.rcParams["text.color"],
                fontweight=mpl.rcParams["axes.titleweight"],
                pad=10,
            )
            right.set_box_aspect(1)
            right.tick_params(axis="y", pad=5)
            right.xaxis.set_major_locator(MaxNLocator(5))
            for label, row in zip(right.get_yticklabels(), chunk.itertuples()):
                if getattr(row, "candidate", False):
                    label.set_fontweight("bold")
    if panel != "contributions" and result.metadata["threshold"] is not None:
        handles, legend_labels = axes[-1].get_legend_handles_labels()
        fig.legend(
            handles,
            legend_labels,
            frameon=False,
            fontsize=mpl.rcParams["legend.fontsize"],
            loc="outside lower center",
            labelcolor=mpl.rcParams["text.color"],
        )
    for axis in axes:
        _tick_fonts(axis)
    if panel in {"both", "contributions"}:
        _labels(left, frame, labels, text)
    if show:
        plt.show()
    return fig


@_with_context
def table(
    result,
    *,
    ax=None,
    show=True,
    decimals=3,
    figsize=None,
    color=None,
    fontsize=None,
    palette="vega20",
    context="paper",
    font_scale=1,
    fontfamily=None,
):
    """Draw a restrained, readable results table; return Axes.

    Args:
        result: EvaluationResult or AuditResult. Values are never normalized.
        ax (Axes, optional): Destination. Cannot accompany figsize.
        show (bool): Display the figure.
        decimals (int): Decimal places for scores, not integer counts.
        figsize (tuple, optional): Figure size in inches.
        color (color): Header accent color.
        fontsize (float, optional): Explicit body text size in points; overrides context.
        context (str or None): Local Seaborn context, or None to inherit active settings.
        font_scale (float): Positive font multiplier relative to context.
        fontfamily (str, optional): Unified family, or inherit the current family.
        palette: Group colors for audit row swatches, as in pl.audit. Ignored
            for evaluation-summary tables.

    Notes:
        Missing values use an em dash. Audit status distinguishes unavailable,
        unscreened, candidate and not-selected groups. Thin black rules on a
        white background separate the header and body without score shading.
    """
    if isinstance(decimals, bool) or not isinstance(decimals, Integral) or decimals < 0:
        raise ValueError("decimals must be a nonnegative integer.")
    color = mpl.rcParams["axes.edgecolor"] if color is None else color
    fontsize = _points(mpl.rcParams["font.size"]) if fontsize is None else fontsize
    if (
        isinstance(fontsize, bool)
        or not isinstance(fontsize, Real)
        or not np.isfinite(fontsize)
        or fontsize <= 0
    ):
        raise ValueError("fontsize must be positive and finite.")
    header_size = _points(mpl.rcParams["xtick.labelsize"])
    to_rgba(color)
    group_colors = {}
    if "itca_plus_gain" in result.table:
        parameter = result.metadata["parameterization"]
        values = result.table.set_index(parameter)[
            ["n_groups", "conditional_itca_sum", "penalty", "itca_plus_gain"]
        ].copy()
        values.columns = [
            "Groups",
            "Conditional\nITCA sum",
            "Penalty",
            "ITCA+\nsplit gain",
        ]
        title, index_label = "Hierarchy selection · nats", (
            "Delta" if parameter == "delta" else "Population ratio"
        )
    elif "node_id" in result.table and "view" in result.table:
        values = result.table.pivot(index="node_id", columns="view", values="itca")
        values.columns = [_view_label(v) + "\nConditional ITCA" for v in values.columns]
        info = result.table.drop_duplicates("node_id").set_index("node_id")
        values["Cells"] = info.n_cells
        values["Status"] = info.reason.replace("", "Evaluated")
        title, index_label = "Within-parent evaluation · nats", "Parent path"
    elif "itca_plus" in result.table:
        source = result.table.copy()
        parameter = result.metadata.get("parameterization", "population_ratio")
        prefix = "δ=" if parameter == "delta" else "p="
        source["label"] = (
            source["groupby"]
            + " · "
            + prefix
            + source[parameter].map(lambda value: f"{value:g}")
        )
        values = source.set_index("label")[
            ["n_groups", "itca", "penalty", "itca_plus", "selected"]
        ]
        values.columns = ["Groups", "ITCA", "Penalty", "ITCA+", "Selected"]
        values["Selected"] = values["Selected"].map({True: "Yes", False: ""})
        title, index_label = "ITCA+ selection · nats", "Grouping · ratio"
    elif "itca" in result.table:
        values = result.table.pivot(index="groupby", columns="view", values="itca")
        views = list(values.columns)
        values.columns = [_view_label(v) + "\nITCA · nats" for v in views]
        sizes = result.table.drop_duplicates("groupby").set_index("groupby")
        values["Groups"] = sizes.n_groups
        values["Cells"] = sizes.n_cells
        title, index_label = "Evaluation summary", "Grouping"
    else:
        source = result.table.set_index("group")
        group_colors = _group_colors(list(source.index), palette)
        values = source[
            ["n_cells", "reference_contribution", "contrast_contribution", "contrast"]
        ].copy()
        values.columns = [
            "Cells",
            _view_label(result.metadata["reference"]) + "\nContribution",
            _view_label(result.metadata["contrast"]) + "\nContribution",
            "Centered\ncontrast",
        ]
        values["Status"] = [
            (
                "Not eligible"
                if not row.eligible
                else (
                    "Candidate"
                    if bool(row.get("candidate", False))
                    else "Not selected" if "candidate" in source else "Not screened"
                )
            )
            for _, row in source.iterrows()
        ]
        title, index_label = "View-specific audit", "Group"
    fig, axes = _axes(
        ax,
        "contrast",
        (
            figsize
            or ((max(8, len(values.columns) * 1.6), max(2.1, len(values) * 0.4 + 1.6)))
            if ax is None
            else figsize
        ),
        len(values),
        max(1, len(values)),
    )
    axis = axes[0]
    nrows, ncols = values.shape
    widths = np.array(
        [1.05]
        + [
            1.35 if "\n" in name or name == "Status" else 0.85
            for name in values.columns
        ]
    )
    edges = np.concatenate([[0], np.cumsum(widths / widths.sum())])
    axis.set(xlim=(0, 1), ylim=(nrows + 1.7, -0.9))
    axis.set_axis_off()
    axis.text(
        0,
        -0.5,
        title,
        fontsize=mpl.rcParams["axes.titlesize"],
        fontweight="bold",
        color=mpl.rcParams["text.color"],
        va="center",
    )
    axis.plot([0, 1], [1, 1], color=color, lw=mpl.rcParams["axes.linewidth"])
    for j, header in enumerate([index_label] + list(values.columns)):
        axis.text(
            edges[j] + 0.016 if j == 0 else (edges[j] + edges[j + 1]) / 2,
            0.45,
            header,
            ha="left" if j == 0 else "center",
            va="center",
            fontsize=header_size,
            color=mpl.rcParams["text.color"],
            linespacing=1.4,
        )
    for i, (index, row) in enumerate(values.iterrows()):
        y = i + 1.5
        if index in group_colors:
            axis.add_patch(
                Rectangle(
                    (0.004, y - 0.19),
                    0.005,
                    0.38,
                    facecolor=group_colors[index],
                    edgecolor="none",
                )
            )
        for j, value in enumerate([index] + list(row)):
            header = index_label if j == 0 else values.columns[j - 1]
            text = (
                "—"
                if pd.isna(value)
                else (
                    f"{int(value):,}"
                    if header in {"Cells", "Groups"}
                    else (
                        f"{value:.{decimals}f}"
                        if isinstance(value, (float, np.floating))
                        else str(value)
                    )
                )
            )
            axis.text(
                edges[j] + 0.016 if j == 0 else (edges[j] + edges[j + 1]) / 2,
                y,
                text,
                ha="left" if j == 0 else "center",
                va="center",
                fontsize=fontsize,
                color=mpl.rcParams["text.color"],
                fontweight="normal",
            )
    axis.plot(
        [0, 1], [nrows + 1, nrows + 1], color=color, lw=mpl.rcParams["axes.linewidth"]
    )
    axis.text(
        0,
        nrows + 1.4,
        (
            "Conditional split objective · not flat partition ITCA+"
            if "itca_plus_gain" in result.table
            else (
                "Within-parent scores · missing values are not zero"
                if "node_id" in result.table and "view" in result.table
                else (
                    "ITCA+ = ITCA − penalty · highest value selected"
                    if "itca_plus" in result.table
                    else (
                        "Raw scores in nats · no composite ranking"
                        if "itca" in result.table
                        else "Contributions and contrasts in nats · descriptive audit"
                    )
                )
            )
        ),
        fontsize=header_size,
        color=mpl.rcParams["text.color"],
        va="center",
    )
    if show:
        plt.show()
    return axis


@_with_context
def hierarchy(
    result,
    *,
    population_ratio=None,
    delta=None,
    palette="vega20",
    ax=None,
    figsize=None,
    show=True,
    context="paper",
    font_scale=1,
    fontfamily=None,
):
    """Plot a selected hierarchy; colors identify its retained terminals.

    Args:
        result (HierarchySelectionResult): Stored pruning, never recomputed.
        population_ratio (float, optional): Choose one stored population ratio.
        delta (float, optional): Alternatively choose one stored delta. One selector
            is required for a multi-parameter result; do not supply both.
        palette: Vega20, mapping from 'level: label' to colors, or None.
        ax (Axes, optional): Existing Matplotlib axis. Cannot accompany figsize.
        figsize (tuple, optional): Figure dimensions in inches.
        show (bool): Display without saving. Returns Axes.
        context, font_scale, fontfamily: Local typography, as in pl.audit.

    Notes:
        Gray edges lead to inactive descendants; bold colored nodes are retained
        terminals. Top-level nodes are fixed roots. No branch length is inferred.
    """
    from ._hierarchy import selected_delta
    from ._results import HierarchySelectionResult

    if not isinstance(result, HierarchySelectionResult):
        raise ValueError(
            "INVALID_RESULT: hierarchy plot requires a hierarchy selection."
        )
    chosen_delta = selected_delta(result, population_ratio, delta)
    nodes = result.nodes
    chosen = result.splits.loc[result.splits.delta == chosen_delta].set_index("node_id")
    children = {node: [] for node in nodes.node_id}
    for row in nodes.itertuples():
        if pd.notna(row.parent_id):
            children[row.parent_id].append(row.node_id)
    roots = nodes.loc[nodes.parent_id.isna(), "node_id"].tolist()
    ordered_leaves, stack = [], roots[::-1]
    while stack:
        node = stack.pop()
        if children[node]:
            stack.extend(children[node][::-1])
        else:
            ordered_leaves.append(node)
    y = {node: i for i, node in enumerate(ordered_leaves)}
    for row in nodes.sort_values("depth", ascending=False, kind="stable").itertuples():
        if children[row.node_id]:
            y[row.node_id] = float(np.mean([y[c] for c in children[row.node_id]]))
    terminals = nodes.loc[nodes.node_id.isin(chosen.index[chosen.terminal])]
    names = sorted(f"{row.level}: {row.label}" for row in terminals.itertuples())
    colors = _group_colors(names, palette)
    if ax is not None and figsize is not None:
        raise ValueError("figsize cannot be used with ax.")
    if ax is None:
        _, ax = plt.subplots(
            figsize=figsize or (6, max(3.2, len(ordered_leaves) * 0.38)),
            layout="constrained",
        )
    _style_axis(ax)
    lookup = nodes.set_index("node_id")
    for row in nodes.itertuples():
        decision = chosen.loc[row.node_id]
        if pd.notna(row.parent_id):
            parent = lookup.loc[row.parent_id]
            ax.plot(
                [parent.depth, row.depth],
                [y[row.parent_id], y[row.node_id]],
                color="black" if decision.active else "#CCCCCC",
                lw=mpl.rcParams["lines.linewidth"],
                zorder=1,
            )
        ax.scatter(
            [row.depth],
            [y[row.node_id]],
            c=[colors[f"{row.level}: {row.label}"] if decision.terminal else "white"],
            edgecolors="black",
            s=36 if decision.terminal else 18,
            zorder=3,
        )
        ax.annotate(
            row.label,
            (row.depth, y[row.node_id]),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center",
            color="black",
            fontsize=mpl.rcParams["xtick.labelsize"],
            fontweight="bold" if decision.terminal else "normal",
        )
    levels = result.metadata["hierarchy"]
    ax.set(
        xticks=range(len(levels)),
        xticklabels=levels,
        yticks=[],
        xlim=(-0.3, len(levels) - 0.7),
        ylim=(len(ordered_leaves) - 0.5, -1),
        title="Selected hierarchy",
    )
    for spine in ax.spines.values():
        spine.set_visible(False)
    _tick_fonts(ax)
    if show:
        plt.show()
    return ax
