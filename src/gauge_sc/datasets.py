"""Synthetic examples and explicit, version-pinned access to real tutorial data."""

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from sklearn.decomposition import PCA

from ._tutorial_data import load_example


def audit_example(*, random_state=0):
    """Generate a controlled two-view simulation with 560 cells and 50 features.

    M0-M3 have distinct molecular programs and indistinguishable QC distributions.
    Q0-Q2 share one molecular program and have distinct QC distributions. Each
    group has 80 cells. This construction illustrates view-specific screening;
    it is not a realistic sequencing simulator or a performance benchmark.

    X contains nonnegative synthetic expression features, not raw counts.
    X_pca is computed before cross-validation. QC-like covariates are generated
    independently, NOT calculated from X; no counts layer is provided.
    obs['simulation_role'] records known design roles for inspection only and
    must not be included as a prediction feature. Returns AnnData; no downloads.
    """
    rng = np.random.default_rng(random_state)
    names = [f"M{i}" for i in range(4)] + [f"Q{i}" for i in range(3)]
    group = np.repeat(np.arange(7), 80)
    program = np.minimum(group, 4)
    profiles = np.ones((5, 50))
    for i in range(5):
        profiles[i, i * 10 : (i + 1) * 10] += 3
    expression = np.maximum(
        profiles[program] + rng.normal(0, 0.6, (len(group), 50)), 0
    ).astype(np.float32)

    # The four M groups share this distribution; Q groups shift independently
    # of molecular expression. Roles are never passed to the evaluator.
    log_depth = rng.normal(9, 0.12, len(group))
    log_detected = rng.normal(7, 0.08, len(group))
    mito = rng.normal(5, 0.5, len(group))
    for i, shift in enumerate([-1.5, 0, 1.5]):
        mask = group == 4 + i
        log_depth[mask] += shift
        log_detected[mask] += shift / 3
        mito[mask] += 10 * (i + 1)
    obs = pd.DataFrame(
        {
            "group": pd.Categorical([names[x] for x in group], categories=names),
            "simulation_role": pd.Categorical(
                np.where(group < 4, "molecular", "qc_associated")
            ),
            "total_counts": np.rint(np.exp(log_depth)),
            "n_genes_by_counts": np.rint(np.exp(log_detected)),
            "pct_counts_mt": np.clip(mito, 0, 100),
        },
        index=[f"synthetic_{i:04d}" for i in range(len(group))],
    )
    adata = AnnData(
        sparse.csr_matrix(expression),
        obs=obs,
        var=pd.DataFrame(index=[f"feature_{i}" for i in range(50)]),
    )
    adata.obsm["X_pca"] = PCA(n_components=15, random_state=random_state).fit_transform(
        expression
    )
    adata.uns["gauge_sc_example"] = {
        "kind": "synthetic",
        "design_version": 3,
        "random_state": random_state,
        "description": "Controlled two-view teaching example; not biological validation.",
        "qc_provenance": "Independent synthetic covariates; not derived from X.",
        "expression_scale": "Nonnegative synthetic features; not counts.",
    }
    return adata
