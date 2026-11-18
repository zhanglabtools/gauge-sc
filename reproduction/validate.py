"""Check numerical relationships in frozen main-figure source data."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from sklearn.metrics import roc_curve

from render import digest, table, verify


def validate_figure4(bundle):
    """Validate both mixtures against the pinned Figure 4 caption and source tables."""
    cells = table(bundle, "fig4b")
    roc = table(bundle, "fig4e")
    null = table(bundle, "fig4f")
    mixtures = {
        "Mixture 1": (
            "fig4c",
            "10x Chromium v2",
            3159,
            (3008, 82, 69),
            (28, 79, 2),
            66,
            6,
            0.967,
        ),
        "Mixture 2": (
            "fig4d",
            "sci-RNA-seq",
            5023,
            (4803, 80, 140),
            (16, 0, 139),
            71,
            5,
            0.903,
        ),
    }
    for frame in (cells, roc, null):
        if set(frame.dataset) != set(mixtures):
            raise ValueError("Figure 4 inputs must cover both manuscript mixtures.")
    if set(roc.record_type) != {"cluster_score", "roc_coordinate"}:
        raise ValueError("Figure 4e requires cluster scores and ROC coordinates.")
    scale = float(cells.contrast_score.abs().quantile(0.99))
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Invalid pooled contrast scale in Figure 4b.")

    summary = {}
    for dataset, (
        panel,
        method,
        n_cells,
        counts,
        mixed,
        n_clusters,
        positives,
        auc_display,
    ) in mixtures.items():
        plotted = cells.loc[cells.dataset == dataset]
        if len(plotted) != n_cells or set(plotted.method) != {method}:
            raise ValueError(f"Figure 4b cell coverage differs for {dataset}.")
        if not np.isfinite(plotted[["umap_1", "umap_2"]]).all().all():
            raise ValueError(
                f"Figure 4b has invalid embedding coordinates for {dataset}."
            )
        if plotted.groupby("cluster").contrast_score.nunique(dropna=False).gt(1).any():
            raise ValueError(
                f"Figure 4b contrast varies within a cluster for {dataset}."
            )

        readout = table(bundle, panel)
        if set(readout.dataset) != {dataset} or set(readout.cell_id) != set(
            plotted.cell_id
        ):
            raise ValueError(f"{panel} cell IDs disagree with Figure 4b.")
        paired = readout.merge(plotted, on="cell_id", validate="one_to_one")
        np.testing.assert_array_equal(paired.cluster_x, paired.cluster_y)
        np.testing.assert_allclose(
            paired.contrast_score_x, paired.contrast_score_y, atol=1e-12, rtol=0
        )
        np.testing.assert_allclose(
            readout.log10_total_counts, np.log10(readout.nUMI + 1), atol=1e-12, rtol=0
        )
        np.testing.assert_allclose(
            readout.log10_detected_genes,
            np.log10(readout.nGene + 1),
            atol=1e-12,
            rtol=0,
        )
        mixed_species = readout.human_fraction_withheld.gt(
            0.15
        ) & readout.human_fraction_withheld.lt(0.85)
        np.testing.assert_array_equal(
            readout.mixed_species_multiplet_withheld, mixed_species
        )
        groups = readout.groupby(
            "candidate_label"
        ).mixed_species_multiplet_withheld.agg(["size", "sum"])
        first = 1 if panel == "fig4c" else 3
        order = ["Other", f"Cand. {first}", f"Cand. {first + 1}"]
        if set(groups.index) != set(order):
            raise ValueError(f"Unexpected candidate labels in {panel}.")
        np.testing.assert_array_equal(groups.loc[order, "size"], counts)
        np.testing.assert_array_equal(groups.loc[order, "sum"], mixed)

        dataset_roc = roc.loc[roc.dataset == dataset]
        clusters = dataset_roc.loc[dataset_roc.record_type == "cluster_score"]
        coordinates = dataset_roc.loc[dataset_roc.record_type == "roc_coordinate"]
        if (
            len(clusters) != n_clusters
            or clusters.duplicated(["method", "cluster"]).any()
            or not clusters.multiplet_enriched.isin([True, False]).all()
            or clusters.multiplet_enriched.sum() != positives
        ):
            raise ValueError(f"Figure 4e cluster coverage differs for {dataset}.")
        fpr, tpr, _ = roc_curve(
            clusters.multiplet_enriched.astype(int), clusters.contrast_score
        )
        np.testing.assert_allclose(
            coordinates.false_positive_rate, fpr, atol=1e-12, rtol=0
        )
        np.testing.assert_allclose(
            coordinates.true_positive_rate, tpr, atol=1e-12, rtol=0
        )
        observed = float(trapezoid(tpr, fpr))
        np.testing.assert_allclose(dataset_roc.auroc, observed, atol=1e-12, rtol=0)
        if round(observed, 3) != auc_display:
            raise ValueError(f"Figure 4e AUROC differs from the caption for {dataset}.")

        draws = null.loc[null.dataset == dataset]
        if (
            len(draws) != 9999
            or set(draws.randomization) != set(range(9999))
            or set(draws.null_type) != {"dual_view"}
            or not draws.auroc.between(0, 1).all()
        ):
            raise ValueError(
                f"Figure 4f matched randomization coverage differs for {dataset}."
            )
        np.testing.assert_allclose(draws.observed_auroc, observed, atol=1e-12, rtol=0)
        exceedances = int(draws.auroc.ge(draws.observed_auroc).sum())
        p_value = (1 + exceedances) / (1 + len(draws))
        np.testing.assert_allclose(
            draws.monte_carlo_p_value, p_value, atol=1e-12, rtol=0
        )
        if p_value != 0.0001:
            raise ValueError(
                f"Figure 4f upper-tail probability differs from the caption for {dataset}."
            )
        summary[dataset] = dict(
            method=method,
            cells=n_cells,
            eligible_cells=int(plotted.contrast_score.notna().sum()),
            groups={
                label: dict(
                    cells=n, mixed_species_cells=k, mixed_species_fraction=k / n
                )
                for label, n, k in zip(order, counts, mixed)
            },
            scored_clusters=n_clusters,
            multiplet_enriched_clusters=positives,
            auroc=observed,
            randomizations=len(draws),
            exceedances=exceedances,
            monte_carlo_p_value=p_value,
        )
    return dict(pooled_absolute_contrast_p99=scale, mixtures=summary)


def validate(bundle, figures=None):
    manifest = verify(bundle)
    figure4 = validate_figure4(bundle)
    checks = [
        "Figure 4b–f: panel inventory and source-file digests verified",
        "Figure 4b–d: cell IDs, cluster contrasts, QC transforms and mixed-species fractions agree",
        "Figure 4e–f: both ROCs and plus-one matched-null probabilities agree with the frozen tables and caption",
    ]
    if figures == [4] or manifest.get("figures") == [4]:
        return dict(
            figures=[4],
            asset_manifest_sha256=digest(bundle / "manifest.json"),
            checks=checks,
            figure4=figure4,
            source_data_checks="passed",
            full_experiment_recomputed=False,
        )
    rare = table(bundle, "fig2i")
    expected_metrics = {"itca_fixed", "itca_adaptive", "f1_macro", "nmi", "silhouette"}
    if set(rare.metric) != expected_metrics or not rare["n"].eq(1).all():
        raise ValueError(
            "Figure 2i must contain five metrics from one deterministic RF schedule per dataset."
        )
    if rare["mean"].isna().any():
        raise ValueError("Figure 2i has unavailable plotted values.")
    consistency = pd.read_csv(bundle / "reference/fig2h_per_dataset.csv")
    for (metric, perturbation), data in consistency.groupby(["metric", "perturbation"]):
        expected = 2 if perturbation == "rare_dropout" else 4
        if data.dataset.nunique() != expected or len(data) != expected:
            raise ValueError(
                f"Figure 2h {metric}/{perturbation} lacks the declared dataset coverage."
            )
    calculated = consistency.groupby(
        ["metric", "perturbation"], as_index=False
    ).consistency.mean()
    recorded = table(bundle, "fig2h")
    merged = calculated.merge(
        recorded,
        on=["metric", "perturbation"],
        suffixes=("_new", "_recorded"),
        validate="one_to_one",
    )
    np.testing.assert_allclose(
        merged.consistency_new, merged.consistency_recorded, atol=1e-12, rtol=0
    )
    checks.append(
        "Figure 2h/i: RF-only single schedules, five metrics, and declared dataset-mean coverage"
    )
    if len(table(bundle, "fig6d_slice")) != 9387:
        raise ValueError("Figure 6d section count differs from Main-v1 (n=9,387).")
    checks.append(
        "Figure 6d: all 9,387 cells in the stated hippocampal section are retained"
    )
    candidates = table(bundle, "fig5d")
    candidates = candidates.loc[candidates.candidate_type == "leiden"].copy()
    path = table(bundle, "fig5f")
    for row in path.itertuples():
        ranked = candidates.assign(
            objective=candidates.mean_aITCA_across_cv_folds
            - row.effective_per_cluster_penalty_delta
            * candidates.number_of_labels_or_clusters
        ).sort_values(
            ["objective", "number_of_labels_or_clusters", "candidate_id"],
            ascending=[False, True, True],
        )
        winner = ranked.iloc[0]
        if winner.candidate_id != row.candidate_id:
            raise ValueError(
                f"Figure 5f winner differs at delta={row.effective_per_cluster_penalty_delta}"
            )
        np.testing.assert_allclose(
            winner.objective, row.selected_candidate_aITCA_MDL_score, atol=1e-12, rtol=0
        )
    checks.append(
        "Figure 5f: all selected candidates and objectives agree with Figure 5d"
    )
    composition = table(bundle, "fig5g")
    np.testing.assert_allclose(
        composition.fraction_of_subclass,
        composition.cell_count / composition.subclass_total_cells,
        atol=1e-12,
        rtol=0,
    )
    checks.append(
        "Figure 5g: displayed fractions retain the complete-subclass denominator"
    )
    for panel in ["fig5h", "fig6j", "fig6k"]:
        data = table(bundle, panel)
        if not data.fraction_expressing.between(-1e-12, 1 + 1e-12).all():
            raise ValueError(f"Invalid marker fractions in {panel}")
    checks.append("Marker dot-plot expression fractions lie in [0, 1]")
    pruning = table(bundle, "fig6b")
    np.testing.assert_array_equal(
        pruning.total_leaves,
        pruning.n_class_leaves + pruning.n_subclass_leaves + pruning.n_supertype_leaves,
    )
    if (
        np.diff(pruning.sort_values("effective_per_cluster_penalty_delta").total_leaves)
        > 0
    ).any():
        raise ValueError("Figure 6b: complexity increases with increasing penalty")
    checks.append(
        "Figure 6b: terminal-level counts add up and pruning complexity is monotone"
    )
    trajectory = table(bundle, "fig7b")
    for _, data in trajectory.groupby("run"):
        data = data.sort_values("trial_index")
        np.testing.assert_allclose(
            data.objective_value.cummax().ffill(),
            data.best_so_far_objective,
            atol=1e-12,
            rtol=0,
        )
    checks.append(
        "Figure 7b: best-so-far trajectories match the recorded trial objectives"
    )
    return dict(
        checks=checks,
        figure4=figure4,
        source_data_checks="passed",
        full_experiment_recomputed=False,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--figures", type=int, nargs=1, choices=[4])
    args = parser.parse_args()
    result = validate(args.bundle, args.figures)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
