"""Local audit validation on a fixed Ding2020 species-blind mixture subset.

Developer-only recipe; requires the validation extra and a local source file.
No data download, species-label selection, distribution, or source modification.
The explicit --output directory receives results and development preview figures.
"""

import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import time

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from scipy import sparse

import gauge_sc as ip


def run(source, output, *, n_cells=2000, seed=0):
    started = time.perf_counter()
    source = Path(source)
    output = Path(output)
    if output.exists():
        raise ValueError("Output already exists; choose a new directory.")
    digest = hashlib.sha256()
    with source.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    backed = ad.read_h5ad(source, backed="r")
    try:
        mask = (backed.obs["Experiment"] == "Mixture1") & (
            backed.obs["Method"] == "10x-Chromium-v2"
        )
        available = np.flatnonzero(mask.to_numpy())
        if not 1 <= n_cells <= len(available):
            raise ValueError("n_cells exceeds the selected sample or is not positive.")
        indices = np.sort(
            np.random.default_rng(seed).choice(available, n_cells, replace=False)
        )
        # Only these whitelisted fields are copied: no withheld species annotations.
        qc = ["nUMI", "nGene", "percent.mito"]
        data = ad.AnnData(
            X=backed.X[indices, :].copy(),
            obs=backed.obs.iloc[indices][["Experiment", "Method", *qc]].copy(),
            var=backed.var.copy(),
        )
        source_processing = str(backed.uns.get("mixture_processing", {}))
        n_available = len(available)
    finally:
        backed.file.close()
    input_genes = data.n_vars
    if (np.asarray(data.X.sum(axis=1)).ravel() <= 0).any():
        raise ValueError(
            "Unexpected zero-count cells; explicit preprocessing revision needed."
        )
    sc.pp.normalize_total(data, target_sum=10000)
    sc.pp.log1p(data)
    sc.pp.highly_variable_genes(data, n_top_genes=2000, flavor="seurat", subset=True)
    sc.pp.scale(data, zero_center=False, max_value=10)
    sc.pp.pca(data, n_comps=30, random_state=seed)
    sc.pp.neighbors(data, n_neighbors=15, use_rep="X_pca", random_state=seed)
    sc.tl.leiden(
        data,
        resolution=1.0,
        random_state=seed,
        flavor="igraph",
        directed=False,
        n_iterations=2,
        key_added="group",
    )
    counts = data.obs["group"].value_counts().sort_index()
    if counts.min() < 5:
        raise ValueError(
            f"Fixed recipe has a group too small for fivefold CV: {counts.to_dict()}; do not silently filter or tune the seed."
        )
    evaluation = ip.evaluate(
        data,
        groupby="group",
        random_state=seed,
        views={
            "molecular": {"obsm": "X_pca"},
            "qc": {"obs": qc, "log1p": ["nUMI", "nGene"]},
        },
    )
    audited = ip.audit(evaluation, reference="molecular", contrast="qc")
    screened = ip.audit(
        evaluation, reference="molecular", contrast="qc", maximum_mirror_ratio=1 / 3
    )
    # Validate the real schema and shared splits without assuming a desired ranking.
    for row in evaluation.table.itertuples():
        pieces = evaluation.groups[evaluation.groups.view == row.view]
        np.testing.assert_allclose(
            pieces.contribution.sum(), row.itca, rtol=0, atol=1e-12
        )
    folds = evaluation.predictions.pivot(
        index="obs_name", columns="view", values="fold"
    )
    np.testing.assert_array_equal(folds.molecular, folds.qc)
    assert not any("withheld" in col for col in data.obs)
    audited.write_to_adata(data)
    sc.tl.umap(data, random_state=seed)
    palette = dict(zip(data.obs["group"].cat.categories, sc.pl.palettes.vega_20_scanpy))
    # Saving is explicit in this validation script; plotting APIs never save.
    output.mkdir(parents=True)
    evaluation.save(output / "evaluation")
    audited.save(output / "audit")
    screened.save(output / "screened")
    data.write_h5ad(output / "local_subset.h5ad", compression="gzip")
    restored = ad.read_h5ad(output / "local_subset.h5ad")
    np.testing.assert_allclose(
        restored.obs["gauge_sc_audit_contrast"],
        data.obs["gauge_sc_audit_contrast"],
        equal_nan=True,
    )
    fig = ip.pl.audit(audited, palette=palette, annotate="all", show=False)
    fig.savefig(output / "audit.png", dpi=160, bbox_inches="tight")
    ax = ip.pl.table(audited, palette=palette, show=False)
    ax.figure.savefig(output / "table.png", dpi=140, bbox_inches="tight")
    fig = sc.pl.umap(
        data,
        color=["group", "gauge_sc_audit_contrast"],
        palette=palette,
        cmap="RdBu_r",
        vcenter=0,
        show=False,
        return_fig=True,
    )
    fig.savefig(output / "umap.png", dpi=140, bbox_inches="tight")
    plt.close("all")
    record = {
        "source_file": str(source.resolve()),
        "source_sha256": digest.hexdigest(),
        "source_processing": source_processing,
        "source_accession": "GSE132044",
        "sample": "Mixture1 / 10x-Chromium-v2",
        "available_cells": n_available,
        "sampled_cells": data.n_obs,
        "sampling": "uniform without replacement, no label stratification",
        "seed": seed,
        "input_genes": input_genes,
        "evaluation_genes": data.n_vars,
        "recipe": "normalize_total 10000; log1p; 2000 Seurat HVGs; scale zero_center=False max_value=10; PCA30; neighbors15; Leiden1.0 igraph 2 iterations; RF100 fivefold pooled ITCA",
        "qc_provenance": "Preserved source nUMI, nGene and percent.mito; not recomputed from HVG subset.",
        "group_counts": {str(k): int(v) for k, v in counts.items()},
        "scores": json.loads(evaluation.table.to_json(orient="records")),
        "eligible_groups": int(audited.table.eligible.sum()),
        "screening_status": screened.metadata["screening_status"],
        "candidate_count": int(screened.table.candidate.sum()),
        "local_subset_bytes": (output / "local_subset.h5ad").stat().st_size,
        "runtime_seconds": round(time.perf_counter() - started, 2),
        "interpretation": "Local real-data integration validation; not paper reproduction or proof that selected groups are artifacts. Withheld species labels were not used. No seed or parameter search.",
        "scanpy_version": version("scanpy"),
    }
    (output / "validation.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cells", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    run(args.input, args.output, n_cells=args.cells, seed=args.seed)
