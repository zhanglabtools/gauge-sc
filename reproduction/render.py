"""Reconstruct Figure 2–7 panels from checksum-verified, frozen source tables."""

import argparse
from html import escape
from importlib.metadata import version
import json
from pathlib import Path
import platform
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np
import pandas as pd
from PIL import Image
import seaborn as sns

from source_asset import COSMX_PANELS, PANEL_COUNTS, SCHEMATICS, digest

METRICS = {
    "itca_fixed": "Ref. ITCA",
    "itca_adaptive": "Cluster ITCA",
    "f1_macro": "Macro-F1",
    "nmi": "NMI",
    "silhouette": "Silhouette",
}
COLORS = dict(zip(METRICS, ["#0072B2", "#56B4E9", "#E69F00", "#CC79A7", "#009E73"]))
FIGURE_SUMMARIES = {
    2: "Perturbation robustness",
    3: "Sequencing-platform comparison",
    4: "Human–mouse mixture audit",
    5: "Allen cortex grouping selection",
    6: "Spatial annotation refinement",
    7: "Workflow search and readouts",
}
PANEL_DESCRIPTIONS = {
    "fig2a": "Unperturbed HIA PBMC UMAP",
    "fig2b": "HIA PBMC after depth downsampling",
    "fig2c": "HIA PBMC after expression mixing",
    "fig2d": "HIA PBMC after cell-count downsampling",
    "fig2e": "Score changes under depth downsampling",
    "fig2f": "Score changes under expression mixing",
    "fig2g": "Score changes under cell-count downsampling",
    "fig2h": "Degradation consistency by metric and perturbation",
    "fig2i": "Score changes under targeted rare-cell-type dropout",
    "fig3b": "QC distributions in PBMC replicate 1",
    "fig3c": "Cell-type UMAPs by sequencing method",
    "fig3d": "Molecular ITCA across effective cluster counts",
    "fig3e": "QC ITCA across effective cluster counts",
    "fig3f": "Molecular and QC ITCA by dataset",
    "fig4b": "Cluster contrast scores in both mixtures",
    "fig4c": "QC and mixed-species fractions in Mixture 1",
    "fig4d": "QC and mixed-species fractions in Mixture 2",
    "fig4e": "ROC curves for multiplet-enriched clusters",
    "fig4f": "AUROCs under matched randomization",
    "fig5a": "Allen cortex UMAP by curated class",
    "fig5b": "IT–ET Glut cells by curated subclass",
    "fig5c": "GAUGE-selected Leiden partition",
    "fig5d": "ITCA across curated and Leiden groupings",
    "fig5e": "Exact and approximate refinement thresholds",
    "fig5f": "Leiden selection across the penalty path",
    "fig5g": "Curated subclass composition of selected clusters",
    "fig5h": "Marker expression in selected IT–ET Glut clusters",
    "fig6b": "Hippocampal terminal labels across penalties",
    "fig6c": "Molecular, spatial and QC ITCA across penalties",
    "fig6d": "Hippocampal section and 3D context",
    "fig6e": "Dentate gyrus map and selected hierarchy",
    "fig6f": "CA1 map and selected hierarchy",
    "fig6g": "CosMx native-label composition by selected grouping",
    "fig6h": "CosMx inhibitory labels and GABA subclasses",
    "fig6i": "CosMx oligodendrocyte and OPC labels",
    "fig6j": "Marker expression in selected GABA subclasses",
    "fig6k": "Oligodendrocyte and OPC marker expression",
    "fig7b": "Best ITCA+ scores across search trials",
    "fig7c": "Grouping sizes across endpoint workflows",
    "fig7d": "Workflow endpoint UMAPs",
    "fig7e": "Workflow endpoints on the larger reference",
}


def verify(bundle):
    manifest = json.loads((bundle / "manifest.json").read_text())
    schema_version = manifest["schema_version"]
    if schema_version not in {1, 2}:
        raise ValueError("Unsupported source asset schema.")
    figures = manifest.get("figures", list(PANEL_COUNTS))
    if figures not in (list(PANEL_COUNTS), [4]):
        raise ValueError("Unsupported figure scope in source asset.")
    expected = {f"fig{f}{chr(97+i)}" for f in figures for i in range(PANEL_COUNTS[f])}
    if schema_version == 2:
        expected -= SCHEMATICS
    ids = [p["id"] for p in manifest["panels"]]
    if len(ids) != len(expected) or set(ids) != expected:
        raise ValueError(f"Panel inventory must match the approved figures {figures}.")
    for panel in manifest["panels"]:
        if "figure" in panel and panel["id"] != f"fig{panel['figure']}{panel['panel']}":
            raise ValueError(f"Inconsistent panel metadata: {panel['id']}")
        if panel["id"].startswith("fig4"):
            tables = [] if panel["id"] == "fig4a" else [f"{panel['id']}.csv"]
            if panel["tables"] != tables:
                raise ValueError(
                    f"Invalid Figure 4 source table mapping: {panel['id']}"
                )
        if panel["id"] not in SCHEMATICS and not panel["tables"]:
            raise ValueError(f"Missing quantitative source: {panel['id']}")
        for name in panel["tables"]:
            if f"source_data/{name}" not in manifest["files"]:
                raise ValueError(f"Source table is not pinned: {name}")
    for name, record in manifest["files"].items():
        path = (bundle / name).resolve()
        if not path.is_relative_to(bundle.resolve()) or not path.is_file():
            raise ValueError(f"Missing or unsafe input: {name}")
        if path.stat().st_size != record["bytes"] or digest(path) != record["sha256"]:
            raise ValueError(f"Checksum mismatch: {name}")
    if schema_version == 1 or manifest.get("panel_numbering") != "paper":
        authority = manifest["authority"]
        if schema_version == 1:
            pinned = manifest["files"].get("reference/" + authority["filename"], {})
            if pinned.get("sha256") != authority["sha256"]:
                raise ValueError("The figure reference checksum changed.")
        elif authority["filename"] != "GAUGE-Main-Figures.pptx" or not re.fullmatch(
            r"[0-9a-f]{64}", authority["sha256"]
        ):
            raise ValueError("Missing figure numbering provenance.")
    return manifest


def table(bundle, stem):
    return pd.read_csv(bundle / "source_data" / f"{stem}.csv")


def cosmx_coordinates(manifest, metadata):
    source = manifest["external_inputs"]["cosmx_metadata"]
    if metadata is None or not metadata.is_file():
        raise ValueError(
            f"Provide --cosmx-metadata {source['filename']} for Figures 6h/i. "
            f"Download instructions: {source['access_url']}"
        )
    if metadata.stat().st_size != source["bytes"] or digest(metadata) != source["sha256"]:
        raise ValueError("The CosMx metadata checksum differs from the frozen input.")
    columns = {"CenterX_global_px": "x", "CenterY_global_px": "y"}
    return pd.read_csv(metadata, usecols=list(columns)).rename(columns=columns)


def axes_grid(count, width=4, height=3):
    columns = min(count, 4)
    rows = int(np.ceil(count / columns))
    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=(columns * width, rows * height),
        squeeze=False,
        layout="constrained",
    )
    for ax in axes.flat[count:]:
        ax.set_visible(False)
    return fig, list(axes.flat[:count])


def categorical_scatter(ax, data, x, y, category, color_column=None, legend=True):
    labels = list(pd.unique(data[category].astype(str)))
    palette = dict(
        zip(labels, sns.color_palette("husl", n_colors=len(labels)).as_hex())
    )
    if color_column:
        palette = dict(
            zip(data[category].astype(str), data[color_column].fillna("#D8D8D8"))
        )
    for label in labels:
        subset = data.loc[data[category].astype(str) == label]
        ax.scatter(
            subset[x],
            subset[y],
            s=2,
            c=palette[label],
            label=label,
            linewidths=0,
            rasterized=True,
        )
    ax.set(xlabel=x, ylabel=y, aspect="equal")
    if legend and len(labels) <= 25:
        ax.legend(fontsize=5, markerscale=3, loc="upper left", bbox_to_anchor=(1, 1))


def heatmap(ax, data, index, columns, values, **kwargs):
    pivot = data.pivot(index=index, columns=columns, values=values)
    sns.heatmap(pivot, ax=ax, **kwargs)


def figure2(panel, data):
    if panel[-1] in "abcd":
        fig, axes = axes_grid(1)
        categorical_scatter(axes[0], data, "UMAP1", "UMAP2", "label")
    elif panel == "fig2h":
        fig, axes = axes_grid(1, 6, 4)
        display = data.assign(metric=data.metric.map(METRICS))
        matrix = display.pivot(
            index="metric", columns="perturbation", values="consistency"
        )
        matrix = matrix.reindex(
            index=list(METRICS.values()),
            columns=["depth", "mix", "cells", "rare_dropout"],
        )
        matrix.columns = [
            "Depth\ndownsampling",
            "Expression\nmixing",
            "Cell-count\ndownsampling",
            "Rare cell\ndropout",
        ]
        sns.heatmap(
            matrix,
            ax=axes[0],
            cmap="RdBu_r",
            vmin=-1,
            vmax=1,
            annot=True,
            fmt=".2f",
            cbar_kws={"label": "Degradation consistency"},
        )
        axes[0].tick_params(axis="y", labelrotation=0)
        axes[0].set_ylabel("")
    else:
        names = list(data.dataset.unique())
        fig, axes = axes_grid(len(names))
        for ax, name in zip(axes, names):
            sub = data.loc[data.dataset == name]
            x = "severity_target" if panel == "fig2i" else "severity"
            if panel == "fig2i":
                sub = sub.loc[sub.control_type == "targeted_rare_drop"]
            for metric, frame in sub.groupby("metric", sort=False):
                frame = frame.sort_values(x)
                ax.plot(
                    frame[x],
                    frame["mean"],
                    "o-",
                    ms=2,
                    color=COLORS[metric],
                    label=METRICS[metric],
                )
                if panel != "fig2i":
                    sd = frame["std"].fillna(0)
                    ax.fill_between(
                        frame[x],
                        frame["mean"] - sd,
                        frame["mean"] + sd,
                        color=COLORS[metric],
                        alpha=0.12,
                    )
            ax.axhline(0, c="grey", ls="--", lw=0.5)
            if panel == "fig2i":
                name = {
                    "HumanImmuneAtlas": "HIA PBMC",
                    "Allen_Cortex_20k": "Allen Cortex",
                }.get(name, name)
                severity = sorted(sub[x].unique())
                ax.set_xticks(severity, [f"{100 * value:g}" for value in severity])
            ax.set(
                title=name,
                xlabel="Target dropout cell mass (%)" if panel == "fig2i" else "Severity",
                ylabel="Score change from baseline",
            )
            ax.legend(fontsize=6)
    return fig


def figure3(panel, data):
    if panel == "fig3b":
        fig, axes = axes_grid(3, 4, 4)
        for ax, col in zip(
            axes,
            [
                "log10_umi_plus_one",
                "log10_detected_genes_plus_one",
                "plotted_mitochondrial_fraction",
            ],
        ):
            sns.violinplot(
                data=data,
                x=col,
                y="display_method",
                ax=ax,
                inner="quart",
                cut=0,
                linewidth=0.5,
            )
    elif panel == "fig3c":
        methods = list(data.sequencing_method.unique())
        fig, axes = axes_grid(len(methods))
        for ax, method in zip(axes, methods):
            categorical_scatter(
                ax,
                data.loc[data.sequencing_method == method],
                "umap_1",
                "umap_2",
                "cell_type",
            )
            ax.set_title(method)
    elif panel in {"fig3d", "fig3e"}:
        fig, axes = axes_grid(1)
        value = "molecular_view_aITCA" if panel == "fig3d" else "qc_view_aITCA"
        for method, sub in data.groupby("sequencing_method", sort=False):
            sub = sub.sort_values("effective_k")
            axes[0].plot(sub.effective_k, sub[value], "o-", ms=3, label=method)
        axes[0].set(
            xlabel="Effective K",
            ylabel="Molecular-view ITCA" if panel == "fig3d" else "QC-view ITCA",
        )
        axes[0].legend(fontsize=6)
    else:
        groups = list(data.dataset_replicate.unique())
        fig, axes = axes_grid(len(groups))
        for ax, group in zip(axes, groups):
            for row in data.loc[data.dataset_replicate == group].itertuples():
                x, y = row.mean_qc_view_aITCA, row.mean_molecular_view_aITCA
                ax.errorbar(
                    x,
                    y,
                    xerr=row.sd_qc_view_aITCA_across_grid,
                    yerr=row.sd_molecular_view_aITCA_across_grid,
                    fmt="o",
                    ms=4,
                    label=row.display_method,
                )
            ax.set(
                title=group,
                xlabel="Mean QC-view ITCA",
                ylabel="Mean molecular-view ITCA",
            )
            ax.legend(fontsize=5)
    return fig


def figure4(panel, data):
    if panel == "fig4b":
        groups = list(data.dataset.unique())
        fig, axes = axes_grid(len(groups))
        bound = float(data.contrast_score.abs().quantile(0.99))
        cmap = LinearSegmentedColormap.from_list(
            "molecular_qc_contrast", ["#0072B2", "#F2F2F2", "#E69F00"]
        ).with_extremes(bad="#E5E8EA")
        norm = TwoSlopeNorm(vmin=-bound, vcenter=0, vmax=bound)
        for ax, name in zip(axes, groups):
            sub = data.loc[data.dataset == name]
            artist = ax.scatter(
                sub.umap_1,
                sub.umap_2,
                c=sub.contrast_score,
                cmap=cmap,
                norm=norm,
                plotnonfinite=True,
                s=3,
                rasterized=True,
            )
            ax.set(
                title=f"{name} · {sub.method.iloc[0]}",
                xlabel="UMAP1",
                ylabel="UMAP2",
                aspect="equal",
            )
        fig.colorbar(artist, ax=axes, label="Contrast score")
    elif panel in {"fig4c", "fig4d"}:
        fig, axes = plt.subplots(2, 2, figsize=(7, 5), layout="constrained")
        candidates = (
            ["Cand. 1", "Cand. 2"] if panel == "fig4c" else ["Cand. 3", "Cand. 4"]
        )
        order = ["Other", *candidates]
        colors = ["#CC79A7", "#009E73"] if panel == "fig4c" else ["#009E73", "#CC79A7"]
        palette = dict(zip(order, ["#D7DCE0", *colors]))
        counts = data.candidate_label.value_counts()
        labels = {
            "log10_total_counts": "log₁₀(UMI + 1)",
            "log10_detected_genes": "log₁₀(genes + 1)",
            "percent.mito": "Mitochondrial counts (%)",
            "mixed_species_multiplet_withheld": "Withheld mixed-species fraction",
        }
        for ax, (col, label) in zip(axes.flat, labels.items()):
            if col == "mixed_species_multiplet_withheld":
                fraction = (
                    data.assign(fraction=data[col].astype(float))
                    .groupby("candidate_label", sort=False)
                    .fraction.mean()
                    .reindex(order)
                )
                bars = ax.bar(order, fraction.values, color=[palette[g] for g in order])
                ax.bar_label(
                    bars, labels=[f"{value:.1%}" for value in fraction], padding=2
                )
                ax.set_ylim(0, 1.12)
            else:
                sns.violinplot(
                    data=data,
                    x="candidate_label",
                    y=col,
                    ax=ax,
                    hue="candidate_label",
                    order=order,
                    hue_order=order,
                    palette=palette,
                    legend=False,
                    inner=None,
                    cut=0,
                    density_norm="width",
                    linewidth=0.6,
                )
                sns.boxplot(
                    data=data,
                    x="candidate_label",
                    y=col,
                    order=order,
                    width=0.20,
                    showfliers=False,
                    color="white",
                    ax=ax,
                )
            ax.set(xlabel="", ylabel=label)
            ax.set_xticks(range(len(order)), [f"{g}\n(n={counts[g]:,})" for g in order])
    elif panel == "fig4e":
        fig, axes = axes_grid(1)
        for name, sub in data.loc[data.record_type == "roc_coordinate"].groupby(
            "dataset", sort=False
        ):
            axes[0].plot(
                sub.false_positive_rate,
                sub.true_positive_rate,
                label=f"{name} (AUROC {sub.auroc.iloc[0]:.3f})",
            )
        axes[0].plot([0, 1], [0, 1], ":", c="grey", lw=0.5)
        axes[0].set(
            xlabel="False-positive rate",
            ylabel="True-positive rate",
            xlim=(0, 1),
            ylim=(0, 1),
        )
        axes[0].legend(fontsize=7)
    else:
        groups = list(data.dataset.unique())
        fig, array = plt.subplots(
            1,
            len(groups),
            figsize=(8, 3),
            sharex=True,
            sharey=True,
            layout="constrained",
        )
        bins = np.linspace(data.auroc.min(), data.auroc.max(), 31)
        axes = np.atleast_1d(array)
        for ax, name in zip(axes, groups):
            sub = data.loc[data.dataset == name]
            ax.hist(sub.auroc, bins=bins, color="#7E9BA5")
            ax.axvline(0.5, c="grey", ls=":", label="Chance")
            ax.axvline(
                sub.observed_auroc.iloc[0], c="#684978", ls="--", label="Observed"
            )
            ax.set(
                title=f"{name} · P = {sub.monte_carlo_p_value.iloc[0]:.4f}",
                xlabel="AUROC under matched randomization",
                ylabel="Randomization count",
                xlim=(0.25, 1.0),
            )
            ax.legend()
    return fig


def figure5(panel, data):
    fig, axes = axes_grid(1, 5, 4)
    ax = axes[0]
    if panel[-1] in "abc":
        category, color = {
            "a": ("class_display", "class_color"),
            "b": ("subclass_display", "subclass_color"),
            "c": ("selected_leiden_cluster", "cluster_color"),
        }[panel[-1]]
        categorical_scatter(ax, data, "umap_1", "umap_2", category, color)
    elif panel == "fig5d":
        for kind, sub in data.groupby("candidate_type", sort=False):
            ax.errorbar(
                sub.number_of_labels_or_clusters,
                sub.mean_aITCA_across_cv_folds,
                yerr=sub.sd_aITCA_across_cv_folds,
                fmt="o",
                ms=3,
                label=kind,
            )
        ax.set(xlabel="Number of labels/clusters (K)", ylabel="ITCA")
        ax.legend()
    elif panel == "fig5e":
        for name, sub in data.groupby("branch_type", sort=False):
            ax.scatter(
                sub.exact_delta_star, sub.approximate_delta_star, s=8, label=name
            )
        maximum = max(data.exact_delta_star.max(), data.approximate_delta_star.max())
        ax.plot([0, maximum], [0, maximum], "--", c="grey")
        ax.set(
            xlabel="Exact transition threshold",
            ylabel="Approximate transition threshold",
        )
        ax.legend()
        inset = ax.inset_axes([0.62, 0.15, 0.33, 0.33])
        inset.scatter(
            data.smaller_fraction_within_parent, data.relative_approximation_error, s=2
        )
        inset.axhline(0.05, c="grey", ls="--")
        inset.set(xlabel="Smaller fraction", ylabel="Relative error")
        inset.tick_params(labelsize=5)
    elif panel == "fig5f":
        ax.plot(data.effective_per_cluster_penalty_delta, data.selected_leiden_k, "o-")
        for level in ["class", "subclass", "supertype"]:
            ax.axhline(
                data[f"{level}_reference_k"].iloc[0], ls="--", lw=0.6, label=level
            )
        ax.set(xlabel="Effective per-cluster penalty", ylabel="Selected Leiden K")
        ax.set_xscale("symlog", linthresh=0.001)
        selected = data.loc[data.is_chosen_main_text_penalty]
        ax.scatter(
            selected.effective_per_cluster_penalty_delta,
            selected.selected_leiden_k,
            s=70,
            facecolors="none",
            edgecolors="#E69F00",
            zorder=5,
        )
        ax.legend()
    elif panel == "fig5g":
        matrix = (
            data.pivot(
                index="subclass_display_order",
                columns="cluster_display_order",
                values="fraction_of_subclass",
            )
            .sort_index()
            .sort_index(axis=1)
        )
        sns.heatmap(matrix, ax=ax, cmap="Blues", vmin=0, vmax=1)
        labels = (
            data.sort_values("subclass_display_order")
            .drop_duplicates("subclass_display_order")
            .subclass_display
        )
        clusters = (
            data.sort_values("cluster_display_order")
            .drop_duplicates("cluster_display_order")
            .selected_leiden_cluster
        )
        ax.set_yticklabels(labels, rotation=0)
        ax.set_xticklabels(clusters, rotation=90)
        ax.set(xlabel="Selected Leiden cluster", ylabel="Curated subclass")
    else:
        artist = ax.scatter(
            data.cluster_display_order,
            data.gene_display_order,
            s=5 + 75 * data.fraction_expressing,
            c=data.plotted_mean_expression_z_score,
            cmap="RdBu_r",
            vmin=-2.5,
            vmax=2.5,
        )
        x = data.sort_values("cluster_display_order").drop_duplicates(
            "cluster_display_order"
        )
        y = data.sort_values("gene_display_order").drop_duplicates("gene_display_order")
        ax.set_xticks(x.cluster_display_order, x.selected_leiden_cluster)
        ax.set_yticks(y.gene_display_order, y.marker_gene)
        ax.set(xlabel="Selected Leiden cluster", ylabel="Marker")
        ax.invert_yaxis()
        fig.colorbar(artist, ax=ax, label="Mean-expression z-score")
    return fig


def figure6(panel, data, bundle):
    if panel in {"fig6b", "fig6c"}:
        fig, axes = axes_grid(1)
        ax = axes[0]
        columns = (
            ["total_leaves"]
            if panel == "fig6b"
            else [
                "molecular_expression_aITCA",
                "spatial_context_aITCA",
                "technical_qc_aITCA",
            ]
        )
        for col in columns:
            ax.plot(
                data.effective_per_cluster_penalty_delta,
                data[col],
                "o-",
                label={
                    "total_leaves": "Total terminal leaves",
                    "molecular_expression_aITCA": "Molecular",
                    "spatial_context_aITCA": "Spatial",
                    "technical_qc_aITCA": "QC",
                }[col],
                color={
                    "molecular_expression_aITCA": "#0072B2",
                    "spatial_context_aITCA": "#8E6BBE",
                    "technical_qc_aITCA": "#E69F00",
                }.get(col),
                ms=3,
            )
        ax.set(
            xlabel="Effective per-cluster penalty",
            ylabel=(
                "Number of terminal leaves" if panel == "fig6b" else "Conditional ITCA"
            ),
        )
        if panel == "fig6b":
            ax.stackplot(
                data.effective_per_cluster_penalty_delta,
                data.n_class_leaves,
                data.n_subclass_leaves,
                data.n_supertype_leaves,
                colors=["#0072B2", "#009E73", "#D55E00"],
                labels=["Class", "Subclass", "Supertype"],
                alpha=0.3,
            )
            for column, color in [
                ("flat_class_k", "#0072B2"),
                ("flat_subclass_k", "#009E73"),
            ]:
                ax.axhline(data[column].iloc[0], color=color, ls="--", lw=0.8)
        selected = data.loc[data.is_selected_delta]
        if not selected.empty:
            ax.axvline(
                selected.effective_per_cluster_penalty_delta.iloc[0],
                color="grey",
                ls="--",
                lw=0.8,
            )
        ax.set_xscale("symlog", linthresh=0.001)
        if panel == "fig6c":
            ax.set(xlim=(-0.00015, 1.05), ylim=(0, 1.05))
        ax.legend(fontsize=6)
    elif panel == "fig6d":
        fig = plt.figure(figsize=(9, 4), layout="constrained")
        ax = fig.add_subplot(121, projection="3d")
        for selected, color in [(False, "#C8C8C8"), (True, "#E69F00")]:
            sub = data.loc[data.selected_slice == selected]
            ax.scatter(
                sub.x_ccf,
                sub.y_ccf,
                sub.z_ccf,
                s=0.4 if not selected else 1,
                c=color,
                alpha=0.25 if not selected else 0.75,
                rasterized=True,
            )
        ax.set(xlabel="x CCF (mm)", ylabel="y CCF (mm)", zlabel="z CCF (mm)")
        ax.set_box_aspect(np.ptp(data[["x_ccf", "y_ccf", "z_ccf"]].to_numpy(), axis=0))
        ax.view_init(elev=22, azim=-53)
        from matplotlib.patches import Rectangle

        section = table(bundle, "fig6d_slice")
        locator = fig.add_subplot(122)
        locator.scatter(
            section.y_ccf, section.z_ccf, s=0.4, c="#888888", rasterized=True
        )
        for name, (y0, y1), (z0, z1), color in [
            ("DG", (2.0, 3.1), (3.9, 7.7), "#8C6D31"),
            ("CA1", (1.4, 2.5), (3.2, 8.1), "#7B4EA3"),
        ]:
            locator.add_patch(
                Rectangle((y0, z0), y1 - y0, z1 - z0, fill=False, ec=color, lw=1)
            )
            locator.text(y0, z0, name, color=color, fontsize=7)
        locator.set(
            xlabel="y CCF (mm)",
            ylabel="z CCF (mm)",
            aspect="equal",
            title=f"Thin section (n={len(section):,})",
        )
        locator.invert_yaxis()

    elif panel in {"fig6e", "fig6f"}:
        fig, axes = axes_grid(3, 4, 4)
        spatial = table(bundle, panel + "_left")
        categorical_scatter(
            axes[0],
            spatial,
            "z_ccf",
            "y_ccf",
            "adaptive_label",
            "identity_color",
            legend=False,
        )
        bounds = {"fig6e": ((7.7, 3.9), (2.0, 3.1)), "fig6f": ((8.1, 3.2), (1.4, 2.5))}
        axes[0].set(xlim=bounds[panel][0], ylim=bounds[panel][1])
        right = table(bundle, panel + "_right")
        if panel == "fig6e":
            levels = {"class": 0, "subclass": 1, "supertype": 2}
            paths = []
            labels = {}
            colors = {}
            for row in right.sort_values("display_order").itertuples():
                leaf = (row.adaptive_level, row.terminal_label)
                path = [leaf]
                labels[leaf] = row.terminal_display
                colors[leaf] = row.identity_color
                if row.adaptive_level != "class":
                    parent_level = (
                        "class" if row.adaptive_level == "subclass" else "subclass"
                    )
                    parent = (parent_level, row.parent_label)
                    path.insert(0, parent)
                    labels[parent] = row.parent_display
                    if parent_level == "subclass":
                        parents = spatial.loc[
                            spatial.subclass == row.parent_label, "class"
                        ].unique()
                        if len(parents) != 1:
                            raise ValueError("Ambiguous spatial taxonomy parent")
                        root = ("class", str(parents[0]))
                        path.insert(0, root)
                        labels[root] = str(parents[0])
                paths.append(path)
            # Group adjacent leaf rows by shared ancestors without changing the taxonomy.
            paths.sort(key=lambda path: tuple(label for _, label in path))
            positions = {}
            descendants = {}
            for i, path in enumerate(paths):
                for node in path:
                    descendants.setdefault(node, []).append(i)
            for node, rows in descendants.items():
                positions[node] = (levels[node[0]], float(np.mean(rows)))
            edges = {(a, b) for path in paths for a, b in zip(path, path[1:])}
            for a, b in sorted(edges):
                x, y = positions[a]
                xx, yy = positions[b]
                axes[1].plot([x, xx], [y, yy], c=colors.get(b, "#999999"), lw=0.8)
            for node, (x, y) in positions.items():
                axes[1].scatter(x, y, s=20, c=colors.get(node, "#BBBBBB"))
                axes[1].text(x + 0.06, y, labels[node], fontsize=5.5, va="center")
            axes[1].set(xlim=(-0.1, 3), ylim=(len(paths) - 0.5, -0.5))
            axes[1].axis("off")
        else:
            positions = {
                row.node_id: (row.x_position, row.y_position)
                for row in right.itertuples()
            }
            for row in right.itertuples():
                x, y = positions[row.node_id]
                if row.parent_node_id in positions:
                    px, py = positions[row.parent_node_id]
                    axes[1].plot([px, x], [py, y], c=row.edge_color, lw=0.8)
                axes[1].scatter(x, y, s=25, c=row.identity_color)
                axes[1].text(x + 0.03, y, row.node_label, fontsize=6)
            axes[1].axis("off")
        pies = table(bundle, "fig6ef_terminal_level_pies")
        pies = pies.loc[pies.panel == panel[-1]]
        axes[2].pie(
            pies.n_active_terminal_leaves,
            labels=pies.adaptive_terminal_level,
            colors=pies.level_color,
            autopct="%.1f%%",
        )
        axes[2].set_title("Active terminal leaves")
    elif panel == "fig6g":
        fig, axes = axes_grid(1, 8, 3)
        matrix = (
            data.pivot(
                index="native_display_order",
                columns="adaptive_terminal_display_order",
                values="fraction_of_native_label",
            )
            .sort_index()
            .sort_index(axis=1)
        )
        sns.heatmap(matrix, ax=axes[0], cmap="Blues", vmin=0, vmax=1)
        rows = data.sort_values("native_display_order").drop_duplicates(
            "native_display_order"
        )
        columns = data.sort_values("adaptive_terminal_display_order").drop_duplicates(
            "adaptive_terminal_display_order"
        )
        axes[0].set_yticklabels(rows.native_label, rotation=0)
        axes[0].set_xticklabels(columns.adaptive_terminal_display, rotation=90)
    elif panel in {"fig6h", "fig6i"}:
        fig, axes = axes_grid(2)
        for ax, side in zip(axes, ["left", "right"]):
            categorical_scatter(
                ax, data, "x", "y", f"{side}_plot_category", f"{side}_plot_color"
            )
            ax.invert_yaxis()
    else:
        fig, axes = axes_grid(1, 6, 4)
        ax = axes[0]
        artist = ax.scatter(
            data.marker_gene_display_order,
            data.group_display_order,
            s=data.dot_size_minimum_points_squared
            + data.dot_size_range_points_squared * data.fraction_expressing,
            c=data.mean_expression_z_score,
            cmap="RdBu_r",
            vmin=-2.5,
            vmax=2.5,
        )
        genes = data.sort_values("marker_gene_display_order").drop_duplicates(
            "marker_gene_display_order"
        )
        groups = data.sort_values("group_display_order").drop_duplicates(
            "group_display_order"
        )
        ax.set_xticks(genes.marker_gene_display_order, genes.marker_gene, rotation=90)
        ax.set_yticks(groups.group_display_order, groups.group)
        ax.invert_yaxis()
        fig.colorbar(artist, ax=ax, label="Mean-expression z-score")
    return fig


def figure7(panel, data):
    if panel == "fig7b":
        fig, axes = axes_grid(1)
        for name, sub in data.groupby("run", sort=False):
            axes[0].plot(sub.trial_index, sub.best_so_far_objective, label=name, lw=0.8)
        axes[0].set(xlabel="Number of BO trials", ylabel="Best-so-far ITCA+")
        axes[0].legend(fontsize=6)
    elif panel == "fig7d":
        names = list(data.workflow.unique())
        fig, axes = axes_grid(len(names))
        for ax, name in zip(axes, names):
            categorical_scatter(
                ax,
                data.loc[data.workflow == name],
                "UMAP1",
                "UMAP2",
                "cluster",
                legend=False,
            )
            ax.set_title(name)
    else:
        fig, axes = axes_grid(3)
        category = "method" if panel == "fig7c" else "endpoint"
        for ax, col in zip(
            axes, ["K", "min_cluster_fraction_pct", "n_clusters_below_0p5pct"]
        ):
            sns.stripplot(data=data, x=category, y=col, ax=ax, jitter=False, size=5)
            ax.tick_params(axis="x", rotation=30)
    return fig


def write_figure_overview(output, records):
    """Combine each source figure's result panels into one overview image."""
    previews = output / "figures"
    previews.mkdir()
    figures = []
    for number, (paper_figure, summary) in enumerate(FIGURE_SUMMARIES.items(), 1):
        panels = [
            record["panel"] for record in records if record["figure"] == paper_figure
        ]
        if not panels:
            continue
        width, height = 640, 420
        columns = min(2, len(panels))
        rows = (len(panels) + columns - 1) // columns
        preview = Image.new("RGB", (columns * width, rows * height), "white")
        for index, panel in enumerate(panels):
            with Image.open(output / f"{panel}.png") as source:
                image = source.convert("RGB")
                image.thumbnail((width - 12, height - 12), Image.LANCZOS)
                row, column = divmod(index, columns)
                preview.paste(
                    image,
                    (
                        column * width + (width - image.width) // 2,
                        row * height + (height - image.height) // 2,
                    ),
                )
        path = previews / f"figure-{number}.png"
        preview.save(path)
        figures.append(
            dict(
                number=number,
                paper_figure=paper_figure,
                summary=summary,
                preview=path.relative_to(output).as_posix(),
                panels=panels,
            )
        )

    cards = "\n".join(
        f'<section><h2>Figure {f["number"]} · {escape(f["summary"])}</h2>'
        f'<a href="{f["preview"]}"><img src="{f["preview"]}" '
        f'alt="{escape(f["summary"])}" loading="lazy"></a></section>'
        for f in figures
    )
    (output / "index.html").write_text(
        '<!doctype html><html lang="en"><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        "<title>GAUGE paper reproduction</title>"
        "<style>body{font:16px system-ui,sans-serif;max-width:76rem;margin:2rem auto;"
        "padding:0 1rem}main{display:grid;gap:1.5rem;"
        "grid-template-columns:repeat(auto-fit,minmax(min(100%,20rem),1fr))}"
        "h2{font-size:1rem}img{display:block;width:100%;max-height:22rem;"
        "object-fit:contain;object-position:top}</style>"
        "<h1>GAUGE paper reproduction</h1><main>" + cards + "</main></html>",
        encoding="utf-8",
    )
    return figures


def main(args):
    bundle = args.bundle.resolve()
    output = args.output.resolve()
    manifest = verify(bundle)
    if args.figures and not set(args.figures) <= set(
        manifest.get("figures", PANEL_COUNTS)
    ):
        raise ValueError("Requested figure is absent from this source asset.")
    coordinates = None
    used_external_inputs = {}
    external_inputs = manifest.get("external_inputs", {})
    if "cosmx_metadata" in external_inputs and (not args.figures or 6 in args.figures):
        coordinates = cosmx_coordinates(manifest, args.cosmx_metadata)
        used_external_inputs["cosmx_metadata"] = external_inputs["cosmx_metadata"]
    if output.exists():
        raise ValueError("Use a new output directory.")
    output.mkdir(parents=True)
    sns.set_theme(style="ticks", font="DejaVu Sans", font_scale=0.75)
    plt.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none"})
    records = []
    for spec in manifest["panels"]:
        panel = spec["id"]
        if panel in SCHEMATICS or (
            args.figures and spec["figure"] not in args.figures
        ):
            continue
        stem = panel + "_left" if panel in {"fig6e", "fig6f"} else panel
        data = table(bundle, stem)
        if data.empty:
            raise ValueError(f"Empty source: {panel}")
        if panel in COSMX_PANELS and coordinates is not None:
            rows = data.pop("metadata_row")
            if (
                not pd.api.types.is_integer_dtype(rows)
                or rows.duplicated().any()
                or not rows.between(0, len(coordinates) - 1).all()
            ):
                raise ValueError(f"Invalid CosMx metadata row mapping in {panel}.")
            data[["x", "y"]] = coordinates.loc[rows, ["x", "y"]].to_numpy()
        renderer = {
            2: figure2,
            3: figure3,
            4: figure4,
            5: figure5,
            6: figure6,
            7: figure7,
        }[spec["figure"]]
        fig = (
            renderer(panel, data, bundle)
            if spec["figure"] == 6
            else renderer(panel, data)
        )
        fig.suptitle(panel, fontsize=11)
        for ext in ["png", "pdf"]:
            fig.savefig(output / f"{panel}.{ext}", dpi=180, bbox_inches="tight")
        plt.close(fig)
        records.append(
            dict(
                panel=panel,
                figure=spec["figure"],
                description=PANEL_DESCRIPTIONS[panel],
                kind=spec["kind"],
                inputs=spec["tables"],
                png_sha256=digest(output / f"{panel}.png"),
            )
        )
        print(f"Rendered {panel}", flush=True)
    figures = write_figure_overview(output, records)
    report = dict(
        mode="source-data reconstruction",
        asset_manifest_sha256=digest(bundle / "manifest.json"),
        external_inputs=used_external_inputs,
        full_experiment_recomputed=False,
        python=platform.python_version(),
        versions={p: version(p) for p in ["numpy", "pandas", "matplotlib", "seaborn"]},
        figures=figures,
        panels=records,
    )
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figures", type=int, nargs="+", choices=range(2, 8))
    parser.add_argument("--cosmx-metadata", type=Path)
    main(parser.parse_args())
