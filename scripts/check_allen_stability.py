"""Compare Allen decisions in fresh processes; never force a particular winner."""

import argparse
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import subprocess
import sys


def worker(args):
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from threadpoolctl import threadpool_info

    import gauge_sc as gauge

    target = args.output
    target.mkdir(parents=True, exist_ok=False)
    adata = gauge.datasets.load_example(
        "allen_cortex", data_dir=args.data, download=False
    )
    sc.pp.pca(adata, n_comps=30, random_state=0)
    sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca", random_state=0)
    keys = []
    for resolution in [0.3, 0.6, 0.9, 1.2, 1.5]:
        key = f"leiden_{resolution}"
        sc.tl.leiden(
            adata,
            resolution=resolution,
            key_added=key,
            random_state=0,
            flavor="igraph",
            directed=False,
            n_iterations=2,
        )
        keys.append(key)
    result = gauge.evaluate(
        adata,
        groupby=keys,
        views={"molecular": {"obsm": "X_pca"}},
        aggregation="fold_mean",
        random_state=args.seed,
    )
    selection = gauge.select(result, view="molecular", delta=[0, 0.01, 0.02, 0.05, 0.1])
    result.save(target / "evaluation")
    selection.save(target / "selection")
    np.savez_compressed(target / "pca.npz", X_pca=adata.obsm["X_pca"])
    adata.obs[keys].to_csv(target / "candidate_labels.csv")
    rows = []
    for delta, candidates in selection.table.groupby("delta", sort=False):
        ranked = candidates.sort_values(
            ["itca_plus", "n_groups", "groupby"], ascending=[False, True, True]
        )
        first, second = ranked.iloc[0], ranked.iloc[1]
        rows.append(
            dict(
                delta=delta,
                winner=first["groupby"],
                n_groups=int(first.n_groups),
                winner_score=float(first.itca_plus),
                runner_up=second["groupby"],
                runner_up_score=float(second.itca_plus),
                score_margin_nats=float(first.itca_plus - second.itca_plus),
            )
        )
    pd.DataFrame(rows).to_csv(target / "decision_margins.csv", index=False)
    record = dict(
        requested_core=os.environ.get("OPENBLAS_CORETYPE", "default"),
        seed=args.seed,
        python=platform.python_version(),
        platform=platform.platform(),
        versions={
            p: version(p)
            for p in ["gauge-sc", "numpy", "scipy", "scanpy", "scikit-learn"]
        },
        threadpools=threadpool_info(),
        pca_sha256=hashlib.sha256(adata.obsm["X_pca"].tobytes()).hexdigest(),
        decisions=rows,
        interpretation="Teaching subset; not the full Figure 5 K=45 experiment.",
    )
    (target / "environment.json").write_text(json.dumps(record, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--coretypes", nargs="+", default=["default", "Haswell"])
    parser.add_argument("--seeds", nargs="+", type=int, default=[0, 1, 2])
    parser.add_argument("--seed", type=int, default=0, help=argparse.SUPPRESS)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.data = args.data.resolve()
    args.output = args.output.resolve()
    if args.worker:
        worker(args)
        return
    if args.output.exists():
        parser.error("Use a new output directory; evidence is never overwritten.")
    args.output.mkdir(parents=True)
    records = []
    for core in args.coretypes:
        if core != "default" and platform.machine().lower() not in {"x86_64", "amd64"}:
            parser.error("Explicit x86 OpenBLAS kernels require a compatible x86 host.")
        for seed in args.seeds:
            env = os.environ.copy()
            env.pop("OPENBLAS_CORETYPE", None)
            if core != "default":
                env["OPENBLAS_CORETYPE"] = core
            env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MPLBACKEND="Agg")
            target = args.output / f"{core}-seed{seed}"
            subprocess.run(
                [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    "--worker",
                    "--data",
                    str(args.data),
                    "--output",
                    str(target),
                    "--seed",
                    str(seed),
                ],
                env=env,
                check=True,
            )
            records.append(json.loads((target / "environment.json").read_text()))
    (args.output / "summary.json").write_text(json.dumps(records, indent=2) + "\n")
    print(f"Recorded {len(records)} fresh-process comparisons in {args.output}")


if __name__ == "__main__":
    main()
