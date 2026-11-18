"""Conditional node evaluation and ancestor-constrained dynamic programming."""

from copy import deepcopy
from numbers import Real
import json

import numpy as np
import pandas as pd

from ._results import HierarchyEvaluationResult, HierarchySelectionResult


def evaluate_hierarchy(
    adata,
    *,
    hierarchy,
    views,
    n_splits,
    random_state,
    aggregation,
    classifier,
    min_child_cells,
):
    from ._evaluate import evaluate, _positive_int, _view, _fingerprint

    if (
        not isinstance(hierarchy, (list, tuple))
        or len(hierarchy) < 2
        or not all(isinstance(x, str) for x in hierarchy)
        or len(set(hierarchy)) != len(hierarchy)
    ):
        raise ValueError(
            "INVALID_HIERARCHY: provide distinct obs columns, coarse to fine."
        )
    _positive_int(n_splits, "n_splits", 2)
    _positive_int(random_state, "random_state", 0)
    _positive_int(min_child_cells, "min_child_cells")
    if aggregation not in {"pooled", "fold_mean"}:
        raise ValueError("INVALID_PARAMETER: aggregation must be pooled or fold_mean.")
    if not adata.obs_names.is_unique or adata.n_obs == 0:
        raise ValueError("INVALID_OBSERVATIONS: unique nonempty observations required.")
    if (
        not isinstance(views, dict)
        or not views
        or not all(isinstance(k, str) and k for k in views)
    ):
        raise ValueError("INVALID_VIEW: supply named view definitions.")
    matrices = {name: _view(adata, name, spec)[0] for name, spec in views.items()}
    labels = pd.DataFrame(index=adata.obs_names)
    for level in hierarchy:
        if level not in adata.obs:
            raise ValueError(f"MISSING_GROUPING: obs[{level!r}] does not exist.")
        raw = adata.obs[level]
        if raw.isna().any():
            raise ValueError(f"MISSING_LABEL: {level!r} has missing assignments.")
        labels[level] = raw.astype(str)
        if labels[level].nunique() != raw.nunique():
            raise ValueError("AMBIGUOUS_LABEL: distinct labels render identically.")
    for parent, child in zip(hierarchy, hierarchy[1:]):
        if (labels.groupby(child)[parent].nunique() > 1).any():
            raise ValueError(
                f"NON_NESTED_HIERARCHY: labels in {child!r} have multiple parents."
            )
    records, positions, children = [], {}, {}
    paths = [[] for _ in range(adata.n_obs)]
    for depth, level in enumerate(hierarchy):
        for i, value in enumerate(labels[level]):
            paths[i] = paths[i] + [value]
        ids = [json.dumps(path, ensure_ascii=False) for path in paths]
        level_positions = {}
        for i, node in enumerate(ids):
            level_positions.setdefault(node, []).append(i)
        for node in sorted(level_positions):
            pos = np.asarray(level_positions[node], dtype=int)
            path = json.loads(node)
            parent = json.dumps(path[:-1], ensure_ascii=False) if depth else None
            positions[node] = pos
            children[node] = []
            if parent is not None:
                children[parent].append(node)
            records.append(
                dict(
                    node_id=node,
                    parent_id=parent,
                    level=level,
                    depth=depth,
                    label=path[-1],
                    n_cells=len(pos),
                )
            )
    nodes = pd.DataFrame(records)
    nodes["n_children"] = nodes.node_id.map(lambda node: len(children[node]))
    tables, predictions, groups, folds, provenance = [], [], [], [], {}
    for node in nodes.itertuples():
        if not children[node.node_id]:
            continue
        pos = positions[node.node_id]
        smallest = min(len(positions[c]) for c in children[node.node_id])
        eligible = node.n_children >= 2 and smallest >= max(n_splits, min_child_cells)
        reason = (
            ""
            if eligible
            else (
                "fewer_than_two_children"
                if node.n_children < 2
                else "insufficient_child_support"
            )
        )
        if eligible:
            local = evaluate(
                adata[pos].copy(),
                groupby=hierarchy[node.depth + 1],
                views=views,
                n_splits=n_splits,
                random_state=random_state,
                aggregation=aggregation,
                classifier=classifier,
            )
            provenance[node.node_id] = local.metadata
            for destination, frame in [
                (predictions, local.predictions),
                (groups, local.groups),
                (folds, local.folds),
            ]:
                destination.append(frame.assign(node_id=node.node_id))
        for view in views:
            tables.append(
                dict(
                    node_id=node.node_id,
                    level=node.level,
                    label=node.label,
                    view=view,
                    n_cells=len(pos),
                    n_children=node.n_children,
                    min_child_size=smallest,
                    eligible=eligible,
                    reason=reason,
                    itca=(
                        float(local.table.loc[local.table.view == view, "itca"].iloc[0])
                        if eligible
                        else np.nan
                    ),
                )
            )
    return HierarchyEvaluationResult(
        table=pd.DataFrame(tables),
        nodes=nodes,
        assignments=pd.DataFrame({"obs_name": list(adata.obs_names), "leaf_id": ids}),
        predictions=(
            pd.concat(predictions, ignore_index=True) if predictions else pd.DataFrame()
        ),
        groups=pd.concat(groups, ignore_index=True) if groups else pd.DataFrame(),
        folds=pd.concat(folds, ignore_index=True) if folds else pd.DataFrame(),
        metadata=dict(
            schema_version=1,
            hierarchy=list(hierarchy),
            views=deepcopy(views),
            view_fingerprints={k: _fingerprint(v) for k, v in matrices.items()},
            n_splits=int(n_splits),
            random_state=int(random_state),
            aggregation=aggregation,
            min_child_cells=int(min_child_cells),
            score_units="nats",
            metric="conditional ITCA",
            scope="Within-parent child prediction; proportions normalized within each parent.",
            root_policy="Top-level groups are fixed roots; no root split is scored.",
            ineligible_policy="Ineligible splits, including unary nodes, block descendants.",
            node_evaluations=provenance,
        ),
    )


def select_hierarchy(result, *, view, ratios, deltas):
    if view not in result.metadata["views"]:
        raise ValueError("MISSING_VIEW: select an evaluated view explicitly.")
    nodes = result.nodes.copy()
    by_id = nodes.set_index("node_id")
    children = {node: [] for node in nodes.node_id}
    roots = []
    for row in nodes.itertuples():
        if pd.isna(row.parent_id):
            roots.append(row.node_id)
        else:
            children[row.parent_id].append(row.node_id)
    scores = result.table.loc[result.table.view == view].set_index("node_id")
    paths, decisions, assignments = [], [], []
    for p, delta in zip(ratios, deltas):
        best, split, immediate = {}, {}, {}
        for row in nodes.sort_values(
            "depth", ascending=False, kind="stable"
        ).itertuples():
            node = row.node_id
            immediate[node] = None
            if children[node] and bool(scores.loc[node, "eligible"]):
                value = float(scores.loc[node, "itca"])
                if not np.isfinite(value) or value < 0:
                    raise ValueError(
                        "INVALID_NODE_SCORE: eligible splits require finite nonnegative ITCA."
                    )
                immediate[node] = value - delta * (len(children[node]) - 1)
                gain = immediate[node] + sum(best[c] for c in children[node])
                best[node], split[node] = max(0.0, gain), gain > 0
            else:
                best[node], split[node] = 0.0, False
        active, terminals = set(), []
        stack = roots[::-1]
        while stack:
            node = stack.pop()
            active.add(node)
            if split[node]:
                stack.extend(children[node][::-1])
            else:
                terminals.append(node)
        accepted = [node for node in active if split[node]]
        total = sum(float(scores.loc[node, "itca"]) for node in accepted)
        penalty = delta * (len(terminals) - len(roots))
        paths.append(
            dict(
                population_ratio=p,
                delta=delta,
                n_groups=len(terminals),
                n_splits=len(accepted),
                conditional_itca_sum=total,
                penalty=penalty,
                itca_plus_gain=total - penalty,
            )
        )
        for row in nodes.itertuples():
            node = row.node_id
            decisions.append(
                dict(
                    population_ratio=p,
                    delta=delta,
                    node_id=node,
                    level=row.level,
                    label=row.label,
                    active=node in active,
                    terminal=node in terminals,
                    accepted=node in accepted,
                    eligible=(
                        bool(scores.loc[node, "eligible"]) if children[node] else False
                    ),
                    reason=(
                        str(scores.loc[node, "reason"]) if children[node] else "leaf"
                    ),
                    immediate_gain=immediate[node],
                    subtree_gain=best[node],
                )
            )
        terminal_set = set(terminals)
        leaf_to_terminal = {}
        for leaf in result.assignments.leaf_id.unique():
            node = leaf
            while node not in terminal_set:
                node = by_id.loc[node, "parent_id"]
            leaf_to_terminal[leaf] = node
        assigned = result.assignments[["obs_name"]].copy()
        assigned["population_ratio"] = p
        assigned["delta"] = delta
        assigned["node_id"] = result.assignments.leaf_id.map(leaf_to_terminal)
        labels = {
            node: f"{by_id.loc[node, 'level']}: {by_id.loc[node, 'label']}"
            for node in terminals
        }
        assigned["label"] = assigned.node_id.map(labels)
        assignments.append(assigned)
    path = pd.DataFrame(paths)
    return HierarchySelectionResult(
        table=path.copy(),
        path=path,
        nodes=nodes,
        splits=pd.DataFrame(decisions),
        assignments=pd.concat(assignments, ignore_index=True),
        original_assignments=result.assignments.copy(),
        metadata=dict(
            schema_version=1,
            criterion="hierarchical ITCA+ split gain",
            view=view,
            population_ratio=ratios,
            delta=deltas,
            parameterization=(
                "delta" if all(p is None for p in ratios) else "population_ratio"
            ),
            hierarchy=list(result.metadata["hierarchy"]),
            formula="sum(conditional_itca - delta * (n_children - 1)) over accepted splits",
            score_units="nats",
            tie_break="Retain a terminal when total subtree gain is zero.",
            interpretation="Optimal ancestor-closed pruning for stored conditional scores; not flat partition ITCA+.",
            evaluation=deepcopy(result.metadata),
        ),
    )


def selected_delta(result, population_ratio=None, delta=None):
    if population_ratio is not None and delta is not None:
        raise ValueError("AMBIGUOUS_PENALTY: choose one parameter.")
    if population_ratio is not None:
        if (
            not isinstance(population_ratio, Real)
            or isinstance(population_ratio, (bool, np.bool_))
            or population_ratio not in result.metadata["population_ratio"]
        ):
            raise ValueError("MISSING_RATIO: choose one evaluated population_ratio.")
        index = result.metadata["population_ratio"].index(population_ratio)
        return result.metadata["delta"][index]
    if delta is None:
        if len(result.metadata["delta"]) != 1:
            raise ValueError(
                "AMBIGUOUS_RATIO: choose an evaluated population_ratio or delta."
            )
        return result.metadata["delta"][0]
    if (
        not isinstance(delta, Real)
        or isinstance(delta, (bool, np.bool_))
        or delta not in result.metadata["delta"]
    ):
        raise ValueError("MISSING_DELTA: choose one evaluated delta.")
    return float(delta)
