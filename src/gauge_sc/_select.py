"""Complexity-calibrated selection from a single fixed-view evaluation."""

from copy import deepcopy
from numbers import Real

import numpy as np
import pandas as pd

from ._results import EvaluationResult, SelectionResult, HierarchyEvaluationResult


def select(result, *, view, population_ratio=None, delta=None):
    """Select fixed candidate groupings using ITCA+ without refitting.

    Args:
        result (EvaluationResult or HierarchyEvaluationResult): Stored flat or node evaluation.
        view (str): Explicit selection view; normally molecular.
        delta (float or sequence[float], optional): Direct nonnegative per-group
            penalty in nats. Mutually exclusive with population_ratio; one is required.
        population_ratio (float or sequence[float], optional): Explicit population-scale
            calibration p in (0, 1/e]. Each group costs delta = -p log(p) nats.
            A sequence evaluates a path in the supplied order. No size filtering.

    Returns:
        SelectionResult or HierarchySelectionResult: Flat table retains every candidate, including
            itca, penalty, itca_plus and selected. path contains selected rows.
            Exact ties favor fewer groups, then lexical grouping-column order.

    Notes:
        ITCA+ = ITCA - delta * n_groups. This is selection among supplied fixed
        candidates. For hierarchy input, dynamic programming maximizes the sum
        of conditional node ITCA minus delta * (children - 1), with ancestor closure.
        Its output is a split-gain objective, not flat partition ITCA+.
        No global optimum or independent biological validation is implied.
        The interval (0, 1/e] makes the penalty monotone in population_ratio.
    """
    if not isinstance(result, (EvaluationResult, HierarchyEvaluationResult)):
        raise ValueError(
            "INVALID_RESULT: select requires a flat or hierarchy EvaluationResult."
        )
    if not isinstance(view, str) or view not in result.metadata["views"]:
        raise ValueError("MISSING_VIEW: select an evaluated view explicitly.")
    if delta is not None and population_ratio is not None:
        raise ValueError(
            "AMBIGUOUS_PENALTY: provide population_ratio or delta, not both."
        )
    parameter = population_ratio if delta is None else delta
    if isinstance(parameter, Real):
        values = [parameter]
    elif isinstance(parameter, (list, tuple, np.ndarray)):
        values = list(parameter) if np.ndim(parameter) == 1 else []
    else:
        values = []
    error = "INVALID_POPULATION_RATIO" if delta is None else "INVALID_DELTA"
    if not values or any(
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, Real)
        or not np.isfinite(value)
        or not (0 < value <= 1 / np.e if delta is None else value >= 0)
        for value in values
    ):
        raise ValueError(
            f"{error}: use population ratios in (0, 1/e] or finite delta >= 0."
        )
    values = [float(value) for value in values]
    if len(set(values)) != len(values):
        raise ValueError(f"{error}: duplicate values are not allowed.")
    ratios = values if delta is None else [None] * len(values)
    deltas = [float(-p * np.log(p)) for p in values] if delta is None else values
    if isinstance(result, HierarchyEvaluationResult):
        from ._hierarchy import select_hierarchy

        return select_hierarchy(result, view=view, ratios=ratios, deltas=deltas)
    frame = result.table.loc[result.table.view == view].copy()
    if frame.empty or frame["groupby"].duplicated().any():
        raise ValueError("INVALID_CANDIDATES: expected one row per candidate/view.")
    numeric = frame[["itca", "n_cells", "n_groups"]].to_numpy(dtype=float)
    if (
        not np.isfinite(numeric).all()
        or (frame.itca < 0).any()
        or (frame.n_groups < 2).any()
        or (frame.n_groups > frame.n_cells).any()
        or (frame.n_groups % 1 != 0).any()
    ):
        raise ValueError("INVALID_CANDIDATES: invalid scores or group counts.")
    cell_set = None
    for row in frame.itertuples():
        predictions = result.predictions.loc[
            (result.predictions["groupby"] == row.groupby)
            & (result.predictions.view == view)
        ]
        cells = set(predictions.obs_name)
        if predictions.obs_name.duplicated().any() or len(cells) != row.n_cells:
            raise ValueError("CELL_MISMATCH: candidate predictions are incomplete.")
        if cell_set is not None and cells != cell_set:
            raise ValueError("CELL_MISMATCH: candidates must use identical cells.")
        cell_set = cells
    tables = []
    for p, delta in zip(ratios, deltas):
        scored = frame.copy()
        scored["population_ratio"] = p
        scored["delta"] = delta
        scored["penalty"] = delta * scored.n_groups
        scored["itca_plus"] = scored.itca - scored.penalty
        ranked = scored.sort_values(
            ["itca_plus", "n_groups", "groupby"],
            ascending=[False, True, True],
            kind="stable",
        )
        winner = ranked.iloc[0]["groupby"]
        scored["selected"] = scored["groupby"] == winner
        tables.append(scored)
    table = pd.concat(tables, ignore_index=True)
    return SelectionResult(
        table=table,
        path=table.loc[table.selected].reset_index(drop=True).copy(),
        metadata=dict(
            schema_version=1,
            criterion="ITCA+",
            formula="itca_plus = itca - delta * n_groups",
            view=view,
            population_ratio=ratios,
            delta=deltas,
            parameterization=(
                "delta" if all(p is None for p in ratios) else "population_ratio"
            ),
            score_units="nats",
            tie_break="Exact score ties: fewer groups, then lexical groupby order.",
            interpretation="Population-scale complexity calibration; not a group-size filter.",
            evaluation=deepcopy(result.metadata),
        ),
    )
