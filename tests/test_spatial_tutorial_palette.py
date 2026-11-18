"""Exercise the actual spatial tutorial cell with unused and numerous categories."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

sc = pytest.importorskip(
    "scanpy", reason="Spatial tutorial requires the tutorial extra"
)


def test_spatial_tutorial_uses_distinct_colors_without_changing_input():
    notebook = json.loads(
        (
            Path(__file__).parents[1] / "notebooks/real_data/12_spatial_hierarchy.ipynb"
        ).read_text()
    )
    cell = next(
        "".join(cell["source"])
        for cell in notebook["cells"]
        if "level_palette =" in "".join(cell["source"])
    )
    adata = AnnData(
        np.ones((25, 1)),
        obs=pd.DataFrame(
            {
                "class": pd.Categorical(
                    [f"c{i % 11}" for i in range(25)],
                    categories=[f"c{i}" for i in range(23)],
                ),
                "adaptive": pd.Categorical([f"a{i}" for i in range(25)]),
                "selected_level": pd.Categorical(["class"] * 25),
                "z_ccf": np.arange(25),
                "y_ccf": np.arange(25),
            },
            index=[f"cell{i}" for i in range(25)],
        ),
    )
    before = adata.copy()
    scope = {"adata": adata, "plt": plt, "sc": sc}
    try:
        exec(compile(cell, "spatial-tutorial-palette", "exec"), scope)
        displayed = scope.get("plot_data", adata)
        for key in ("class", "adaptive"):
            observed = set(displayed.obs[key].astype(str))
            colors = dict(
                zip(displayed.obs[key].cat.categories, displayed.uns[f"{key}_colors"])
            )
            assert len({colors[label] for label in observed}) == len(observed)
        pd.testing.assert_frame_equal(adata.obs, before.obs)
        assert set(adata.obsm) == set(before.obsm)
        assert set(adata.uns) == set(before.uns)
    finally:
        plt.close("all")
