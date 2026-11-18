"""Prepare fixed local manuscript-data subsets; never download or publish data."""

import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def subset(source, rows, columns):
    # Reopen backed inputs for explicit row reads; never copy all expression data.
    a = ad.read_h5ad(source, backed="r")
    try:
        return ad.AnnData(
            X=a.X[rows, :].copy(),
            obs=a.obs.iloc[rows][columns].copy(),
            var=a.var.copy(),
        )
    finally:
        a.file.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--research-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--allen-input",
        type=Path,
        help="Override the Allen source path, e.g. a local NAS copy.",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    records = {}
    specifications = [
        ("ding_pbmc", "Ding2020_PBMC/Ding2020_PBMC1.h5ad", "GSE132044"),
        (
            "ding_mixture",
            "Ding2020_Mixture/Ding2020_Mixture_processed.h5ad",
            "GSE132044",
        ),
        (
            "abc_hpf",
            "abc_hpf_pilot/hpf_full_log2.h5ad",
            "Allen ABC MERFISH; source metadata retained locally",
        ),
        (
            "vasculature",
            "exp05_human_brain_vasculature/processed/waelchli_brain_sorted_endothelial_cellvoyager_masked_50k_seed0.h5ad",
            "GSE256493",
        ),
    ]
    specifications.append(
        (
            "allen_cortex",
            "Allen_Cortex_20k/allen_cortex_20k.h5ad",
            "Allen whole mouse brain isocortex reference",
        )
    )
    for name, relative, accession in specifications:
        rng = np.random.default_rng(0)
        source = args.research_root / "data" / relative
        if name == "allen_cortex" and args.allen_input is not None:
            source = args.allen_input
        a = ad.read_h5ad(source, backed="r")
        obs = a.obs.copy()
        shape = a.shape
        a.file.close()
        if name == "ding_pbmc":
            parts = []
            for method in ["10x Chromium (v3)", "Drop-seq"]:
                pool = np.flatnonzero(
                    ((obs.Method == method) & (obs.CellType != "Unassigned")).to_numpy()
                )
                if len(pool) < 1000:
                    raise ValueError(f"Too few cells for prespecified method: {method}")
                parts.extend(rng.choice(pool, 1000, replace=False))
            rows = np.sort(parts)
            columns = ["Method", "CellType", "nUMI", "nGene", "percent.mito"]
            policy = "1000 cells per prespecified method, uniform without replacement; no cell-type balancing"
            scale = "source expression used as supplied; see upstream preprocessing"
        elif name == "ding_mixture":
            rows = np.flatnonzero(
                (
                    (obs.Experiment == "Mixture1") & (obs.Method == "10x-Chromium-v2")
                ).to_numpy()
            )
            columns = [
                "nUMI",
                "nGene",
                "percent.mito",
                "mixed_species_multiplet_withheld",
            ]
            policy = "all cells in prespecified Mixture1 / 10x Chromium v2; species labels never used for sampling"
            scale = "species-blind shared-gene counts; source QC retained"
        elif name == "abc_hpf":
            pool = np.flatnonzero(
                (
                    (obs.x_ccf >= 6.8)
                    & (obs.x_ccf <= 7.0)
                    & obs.parcellation_structure.isin(["DG", "CA1", "CA2", "CA3"])
                ).to_numpy()
            )
            rows = np.sort(rng.choice(pool, min(3000, len(pool)), replace=False))
            columns = [
                "class",
                "subclass",
                "supertype",
                "x_ccf",
                "y_ccf",
                "z_ccf",
                "parcellation_structure",
            ]
            policy = "uniform sample of up to 3000 cells in prespecified HPF slice x_CCF 6.8–7.0 mm, DG/CA1/CA2/CA3; no taxonomy balancing"
            scale = (
                "provided log2 expression over 500 genes; no additional log transform"
            )
        elif name == "allen_cortex":
            rows = np.sort(rng.choice(shape[0], 3000, replace=False))
            columns = ["class", "subclass", "supertype"]
            policy = "3000 cells uniformly sampled without replacement from the fixed 20000-cell isocortical reference; no label balancing"
            scale = "normalize total counts to 10000 over all source genes, log1p, then subset to 2000 Seurat HVGs; no second normalization or log transform"
        else:
            rows = np.sort(rng.choice(shape[0], 2500, replace=False))
            columns = []
            policy = "2500 cells uniformly sampled from fixed masked 50000-cell search input; donor and biological labels excluded"
            scale = "counts"
        data = subset(source, rows, columns)
        if name in {"vasculature", "allen_cortex"}:
            data.obs["total_counts"] = np.asarray(data.X.sum(axis=1)).ravel()
            data.obs["n_genes_by_counts"] = np.asarray((data.X > 0).sum(axis=1)).ravel()
        if name == "allen_cortex":
            import scanpy as sc

            sc.pp.normalize_total(data, target_sum=10000)
            sc.pp.log1p(data)
            sc.pp.highly_variable_genes(
                data, n_top_genes=2000, flavor="seurat", subset=True
            )
        data.uns["tutorial_provenance"] = dict(
            source_accession=accession, sampling=policy, seed=0, expression_scale=scale
        )
        destination = args.output / f"{name}.h5ad"
        data.write_h5ad(destination, compression="gzip")
        records[name] = dict(
            filename=destination.name,
            sha256=sha256(destination),
            source=relative,
            source_sha256=sha256(source),
            accession=accession,
            n_obs=data.n_obs,
            n_vars=data.n_vars,
            bytes=destination.stat().st_size,
            sampling=policy,
            expression_scale=scale,
        )
        print(name, data.shape, destination.stat().st_size, flush=True)
    (args.output / "manifest.json").write_text(
        json.dumps(dict(version=1, seed=0, datasets=records), indent=2)
    )


if __name__ == "__main__":
    main()
