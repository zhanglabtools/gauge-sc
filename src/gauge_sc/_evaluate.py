"""Within-dataset, fixed-assignment cross-validation."""

import hashlib
from importlib.metadata import version
import json
from numbers import Integral

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from ._results import EvaluationResult


def _positive_int(value, name, minimum=1):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f"INVALID_PARAMETER: {name} must be an integer >= {minimum}.")


def _view(adata, name, spec):
    if not isinstance(spec, dict):
        raise ValueError(f"INVALID_VIEW: {name!r} must be a dictionary.")
    sources = set(spec) & {"obsm", "obs", "layer", "X"}
    if len(sources) != 1 or set(spec) - sources - {"log1p"}:
        raise ValueError(
            f"INVALID_VIEW: {name!r} needs exactly one of obsm, obs, layer, X; optional log1p applies to obs columns."
        )
    source = next(iter(sources))
    key = spec[source]
    if source == "obs":
        if (
            not isinstance(key, list)
            or not key
            or not all(isinstance(k, str) for k in key)
            or len(set(key)) != len(key)
        ):
            raise ValueError(
                "INVALID_VIEW: obs must be a nonempty list of distinct numeric columns."
            )
        missing = set(key) - set(adata.obs.columns)
        if missing:
            raise ValueError(f"MISSING_COLUMN: {sorted(missing)} in view {name!r}.")
        try:
            X = adata.obs[key].to_numpy(dtype=float, copy=True)
        except (ValueError, TypeError) as exc:
            raise ValueError(
                "NONNUMERIC_VIEW: obs features must be numeric in this release."
            ) from exc
        logs = spec.get("log1p", [])
        if (
            not isinstance(logs, list)
            or set(logs) - set(key)
            or len(logs) != len(set(logs))
        ):
            raise ValueError(
                "INVALID_TRANSFORM: log1p must list distinct selected obs columns."
            )
        for column in logs:
            i = key.index(column)
            if np.any(X[:, i] < 0):
                raise ValueError(
                    f"INVALID_TRANSFORM: {column!r} contains negative values."
                )
            X[:, i] = np.log1p(X[:, i])
    else:
        if "log1p" in spec:
            raise ValueError(
                "INVALID_TRANSFORM: log1p configuration applies only to obs views."
            )
        if source == "X":
            if key is not True:
                raise ValueError(
                    "INVALID_VIEW: use {'X': True} to explicitly select adata.X."
                )
            X = adata.X
        else:
            container = adata.obsm if source == "obsm" else adata.layers
            if not isinstance(key, str) or key not in container:
                raise ValueError(f"MISSING_VIEW: {source}[{key!r}] is unavailable.")
            X = container[key]
        try:
            X = (
                X.astype(np.float64).tocsr()
                if sparse.issparse(X)
                else np.asarray(X, dtype=float)
            )
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "NONNUMERIC_VIEW: features must be a numeric matrix."
            ) from exc
    if X.ndim != 2 or X.shape[0] != adata.n_obs or X.shape[1] == 0:
        raise ValueError("INVALID_SHAPE: view must have shape (n_obs, n_features>0).")
    values = X.data if sparse.issparse(X) else X
    if not np.isfinite(values).all():
        raise ValueError(
            f"NONFINITE_VIEW: {name!r} contains missing or infinite values; preprocess explicitly."
        )
    return X, source == "obs"


def _fingerprint(X):
    """Hash supplied feature values without densifying sparse matrices."""
    digest = hashlib.sha256(str(X.shape).encode())
    arrays = (X.indptr, X.indices, X.data) if sparse.issparse(X) else (X,)
    for array in arrays:
        array = np.ascontiguousarray(array)
        digest.update(str(array.dtype).encode())
        digest.update(memoryview(array).cast("B"))
    return digest.hexdigest()


def evaluate(
    adata,
    *,
    groupby=None,
    views,
    n_splits=5,
    random_state=0,
    aggregation="pooled",
    classifier=None,
    hierarchy=None,
    min_child_cells=10,
):
    """Evaluate fixed groupings from explicitly selected data views.

    Args:
        adata (AnnData): Cells by features; never mutated. No automatic filtering.
        groupby (str or list[str], optional): Existing assignment columns in obs.
            Required for flat evaluation; mutually exclusive with hierarchy.
        hierarchy (list[str], optional): Nested obs columns, coarse to fine.
            Returns conditional within-parent scores instead of flat scores.
        min_child_cells (int): Hierarchy only; minimum cells in every child (10).
            A split also needs n_splits cells per child and at least two children.
            Ineligible nodes block descendant selection; no observations are dropped.
        views (dict): Named source declarations: {'obsm': key}, {'obs': [columns]},
            {'layer': key}, or {'X': True}. Numeric obs views are standardized
            within training folds; optional log1p names nonnegative obs columns.
        n_splits (int): Stratified fold count, default 5. Each group needs at
            least this many observations; smaller groups cause an explicit error.
        random_state (int): Nonnegative split and default-classifier seed.
        aggregation (str): 'pooled' OOF recall (audit default), or 'fold_mean'.
        classifier (estimator, optional): Cloneable sklearn classifier. Default:
            100-tree random forest, one job. Set custom
            estimator random states explicitly for reproducibility.

    Returns:
        EvaluationResult or HierarchyEvaluationResult: Summary, group, fold and
            OOF prediction tables. Hierarchy results also retain nodes/assignments;
            their proportions are normalized within parents, not globally.
            Flat group proportions are fixed over all input observations. Different
            views share folds within each grouping. Distinct groupings have
            their own stratified folds. Scores are in nats.

    Notes:
        This evaluates fixed labels, not independently regenerated clusters.
        Supplied embeddings are treated as fixed; they are not refitted in folds.
    """
    if hierarchy is not None:
        if groupby is not None:
            raise ValueError(
                "AMBIGUOUS_GROUPING: choose groupby or hierarchy, not both."
            )
        from ._hierarchy import evaluate_hierarchy

        return evaluate_hierarchy(
            adata,
            hierarchy=hierarchy,
            views=views,
            n_splits=n_splits,
            random_state=random_state,
            aggregation=aggregation,
            classifier=classifier,
            min_child_cells=min_child_cells,
        )
    _positive_int(n_splits, "n_splits", 2)
    _positive_int(random_state, "random_state", 0)
    if aggregation not in {"pooled", "fold_mean"}:
        raise ValueError("INVALID_AGGREGATION: choose pooled or fold_mean.")
    if not adata.obs_names.is_unique:
        raise ValueError("DUPLICATE_CELLS: obs_names must be unique.")
    names = [groupby] if isinstance(groupby, str) else groupby
    if (
        not isinstance(names, list)
        or not names
        or not all(isinstance(k, str) for k in names)
        or len(set(names)) != len(names)
    ):
        raise ValueError(
            "INVALID_GROUPING: groupby must be a string or list of distinct strings."
        )
    if (
        not isinstance(views, dict)
        or not views
        or not all(isinstance(k, str) and k for k in views)
    ):
        raise ValueError("INVALID_VIEWS: provide named view declarations.")
    matrices = {name: _view(adata, name, spec) for name, spec in views.items()}
    base = (
        classifier
        if classifier is not None
        else RandomForestClassifier(
            n_estimators=100, random_state=random_state, n_jobs=1
        )
    )
    prepared = {}
    # Validate all candidates before any expensive fitting.
    for name in names:
        if name not in adata.obs:
            raise ValueError(f"MISSING_GROUPING: obs[{name!r}] does not exist.")
        labels = adata.obs[name]
        if labels.isna().any():
            raise ValueError(f"MISSING_LABEL: {name!r} contains missing assignments.")
        rendered = labels.astype(str)
        if rendered.nunique() != labels.nunique():
            raise ValueError(
                "AMBIGUOUS_LABEL: distinct labels share the same string representation."
            )
        y, categories = pd.factorize(rendered, sort=True)
        counts = np.bincount(y)
        if len(counts) < 2:
            raise ValueError(
                "SINGLE_GROUP: evaluation requires at least two nonempty groups."
            )
        if counts.min() < n_splits:
            small = categories[counts < n_splits].tolist()
            raise ValueError(
                f"INSUFFICIENT_GROUP_SIZE: {name!r} groups {small} have fewer than {n_splits} cells; lower n_splits or explicitly revise the input grouping."
            )
        prepared[name] = (y, categories, counts)
    summaries, groups, folds, predictions = [], [], [], []
    for name, (y, categories, counts) in prepared.items():
        splits = list(
            StratifiedKFold(
                n_splits=n_splits, shuffle=True, random_state=random_state
            ).split(np.zeros(len(y)), y)
        )
        p = counts / len(y)
        mass = -p * np.log(p)
        for view, (X, standardize) in matrices.items():
            predicted, fold_id = np.empty(len(y), dtype=int), np.empty(
                len(y), dtype=int
            )
            recalls, supports = [], []
            for fold, (train, test) in enumerate(splits):
                estimator = clone(base)
                if standardize:
                    estimator = make_pipeline(StandardScaler(), estimator)
                estimator.fit(X[train], y[train])
                pred = estimator.predict(X[test])
                if not np.isin(pred, np.arange(len(categories))).all():
                    raise ValueError(
                        "INVALID_PREDICTION: classifier returned unknown group labels."
                    )
                predicted[test], fold_id[test] = pred, fold
                support = np.bincount(y[test], minlength=len(categories))
                correct = np.bincount(
                    y[test][pred == y[test]], minlength=len(categories)
                )
                recall = correct / support
                recalls.append(recall)
                supports.append(support)
                for i, label in enumerate(categories):
                    folds.append(
                        dict(
                            groupby=name,
                            view=view,
                            group=label,
                            fold=fold,
                            n_test=int(support[i]),
                            recall=float(recall[i]),
                        )
                    )
            pooled = np.bincount(y[predicted == y], minlength=len(categories)) / counts
            mean = np.mean(recalls, axis=0)
            chosen = pooled if aggregation == "pooled" else mean
            for i, label in enumerate(categories):
                groups.append(
                    dict(
                        groupby=name,
                        view=view,
                        group=label,
                        n_cells=int(counts[i]),
                        proportion=float(p[i]),
                        information_mass=float(mass[i]),
                        recall=float(chosen[i]),
                        recall_pooled=float(pooled[i]),
                        recall_fold_mean=float(mean[i]),
                        contribution=float(mass[i] * chosen[i]),
                        min_test_cells=int(np.min(supports, axis=0)[i]),
                    )
                )
            summaries.append(
                dict(
                    groupby=name,
                    view=view,
                    n_cells=len(y),
                    n_groups=len(categories),
                    entropy=float(mass.sum()),
                    effective_groups=float(np.exp(mass.sum())),
                    itca=float(mass @ chosen),
                    fold_score_sd=float(np.std(np.asarray(recalls) @ mass, ddof=1)),
                )
            )
            predictions.extend(
                dict(
                    groupby=name,
                    view=view,
                    obs_name=obs,
                    group=categories[true],
                    predicted_group=categories[pred],
                    fold=int(fold),
                )
                for obs, true, pred, fold in zip(adata.obs_names, y, predicted, fold_id)
            )
    metadata = dict(
        schema_version=1,
        package_version="0.1.0a2",
        aggregation=aggregation,
        n_splits=int(n_splits),
        random_state=int(random_state),
        views=json.loads(json.dumps(views)),
        view_fingerprints={
            name: _fingerprint(matrix) for name, (matrix, _) in matrices.items()
        },
        grouping_fingerprints={
            name: hashlib.sha256(
                json.dumps([str(categories[i]) for i in y]).encode()
            ).hexdigest()
            for name, (y, categories, _) in prepared.items()
        },
        classifier=repr(base),
        classifier_class=f"{type(base).__module__}.{type(base).__name__}",
        versions={
            pkg: version(pkg) for pkg in ["numpy", "pandas", "scikit-learn", "anndata"]
        },
        obs_names_sha256=hashlib.sha256(
            json.dumps(list(adata.obs_names)).encode()
        ).hexdigest(),
        metric="ITCA",
        score_units="nats",
        fold_score_sd_description="Variation across overlapping training folds; not a confidence interval",
        evaluation_scope="Fixed assignments and fixed input representations; within-dataset prediction",
    )
    return EvaluationResult(
        pd.DataFrame(summaries),
        pd.DataFrame(groups),
        pd.DataFrame(folds),
        pd.DataFrame(predictions),
        metadata,
    )
