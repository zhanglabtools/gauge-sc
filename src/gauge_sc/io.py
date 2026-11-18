"""Load plain JSON results for inspection, plotting and further computation."""

import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from ._results import (
    AuditResult,
    EvaluationResult,
    HierarchyEvaluationResult,
    HierarchySelectionResult,
    SelectionResult,
)

_RESULT_TYPES = {
    "evaluation": (EvaluationResult, ("table", "groups", "folds", "predictions")),
    "audit": (AuditResult, ("table", "assignments")),
    "selection": (SelectionResult, ("table", "path")),
    "hierarchy_evaluation": (
        HierarchyEvaluationResult,
        ("table", "nodes", "assignments", "predictions", "groups", "folds"),
    ),
    "hierarchy_selection": (
        HierarchySelectionResult,
        ("table", "path", "nodes", "splits", "assignments", "original_assignments"),
    ),
}
_REQUIRED = {
    "evaluation": {
        "table": {"groupby", "view", "n_cells", "n_groups", "itca"},
        "groups": {
            "groupby",
            "view",
            "group",
            "n_cells",
            "proportion",
            "information_mass",
            "contribution",
            "recall",
            "min_test_cells",
        },
        "folds": {"groupby", "view", "group", "fold", "n_test", "recall"},
        "predictions": {
            "groupby",
            "view",
            "group",
            "obs_name",
            "predicted_group",
            "fold",
        },
    },
    "audit": {
        "table": {
            "group",
            "contrast",
            "difference",
            "eligible",
            "reason",
            "reference_contribution",
            "contrast_contribution",
        },
        "assignments": {"obs_name", "group"},
    },
    "selection": {
        "table": {
            "groupby",
            "view",
            "itca",
            "delta",
            "penalty",
            "itca_plus",
            "selected",
        },
        "path": {"groupby", "delta", "itca_plus"},
    },
    "hierarchy_evaluation": {
        "table": {"node_id", "view", "n_children", "itca", "eligible", "reason"},
        "nodes": {"node_id", "parent_id", "level", "label", "n_cells", "n_children"},
        "assignments": {"obs_name", "leaf_id"},
    },
    "hierarchy_selection": {
        "table": {"delta", "n_groups", "itca_plus_gain"},
        "path": {"delta", "n_groups", "itca_plus_gain"},
        "nodes": {"node_id", "parent_id", "level", "label", "n_cells", "n_children"},
        "splits": {
            "node_id",
            "delta",
            "active",
            "terminal",
            "accepted",
            "eligible",
            "reason",
        },
        "assignments": {"obs_name", "delta", "node_id", "label"},
        "original_assignments": {"obs_name", "leaf_id"},
    },
}
_METADATA_KEYS = {
    "evaluation": {"views", "n_splits", "aggregation"},
    "hierarchy_evaluation": {"views", "hierarchy", "n_splits", "aggregation"},
    "audit": {"groupby", "reference", "contrast", "maximum_mirror_ratio", "evaluation"},
    "selection": {"view", "delta", "population_ratio", "evaluation"},
    "hierarchy_selection": {
        "view",
        "hierarchy",
        "delta",
        "population_ratio",
        "evaluation",
    },
}

_IDENTITIES = {
    "obs_name",
    "group",
    "predicted_group",
    "groupby",
    "view",
    "node_id",
    "parent_id",
    "leaf_id",
    "level",
    "label",
    "reason",
}
_FLAGS = {"eligible", "candidate", "selected", "active", "terminal", "accepted"}
_STR_RETAINS_NULLS = bool(pd.Series([None], dtype=object).astype("str").isna().iloc[0])
_DTYPE = re.compile(r"(?:object|str|string|bool|u?int(?:8|16|32|64)|float(?:32|64))\Z")


def _invalid(message):
    raise ValueError(f"INVALID_RESULT_FILE: {message}")


def _nonfinite(value):
    _invalid(f"nonfinite JSON constant {value!r}; undefined values must be null.")


def _frame(records, schema, name):
    if not isinstance(records, list) or any(not isinstance(r, dict) for r in records):
        _invalid(f"table {name!r} must contain JSON records.")
    if records and any(set(r) != set(records[0]) for r in records):
        _invalid(f"table {name!r} has inconsistent columns.")
    if schema is not None:
        if not isinstance(schema, dict) or any(
            not isinstance(k, str) or not isinstance(v, str) or not _DTYPE.fullmatch(v)
            for k, v in schema.items()
        ):
            _invalid(f"table {name!r} has an unsupported column schema.")
        if records and set(records[0]) != set(schema):
            _invalid(f"table {name!r} does not match its column schema.")
    frame = pd.DataFrame(records, columns=list(schema) if schema else None)
    for column in frame:
        values = frame[column].dropna().tolist()
        if column in _IDENTITIES and any(not isinstance(v, str) for v in values):
            _invalid(f"{name}.{column} must retain string identities.")
        if column in _FLAGS and (
            frame[column].isna().any() or any(type(v) is not bool for v in values)
        ):
            _invalid(f"{name}.{column} must contain booleans.")
        if schema is None:
            continue
        dtype = schema[column]
        if dtype == "bool" and (
            frame[column].isna().any() or any(type(v) is not bool for v in values)
        ):
            _invalid(f"{name}.{column} must contain booleans.")
        if dtype.startswith(("int", "uint", "float")):
            if any(isinstance(v, (str, bool)) for v in values):
                _invalid(f"{name}.{column} must contain numbers.")
            if dtype.startswith(("int", "uint")):
                numeric = frame[column].to_numpy(dtype=float)
                if not np.isfinite(numeric).all() or (numeric % 1 != 0).any():
                    _invalid(f"{name}.{column} must contain finite integers.")
                bounds = np.iinfo(dtype)
                if any(v < bounds.min or v > bounds.max for v in values):
                    _invalid(f"{name}.{column} exceeds the range of {dtype}.")
        if dtype in ("str", "string"):
            if any(not isinstance(v, str) for v in values):
                _invalid(f"{name}.{column} must contain strings.")
            # pandas 3's str dtype retains missing values; pandas 1/2's astype(str)
            # turns null parents into literal strings and corrupts hierarchy roots.
            if dtype == "str" and not _STR_RETAINS_NULLS:
                dtype = "object"
        try:
            frame[column] = frame[column].astype(dtype)
        except (ValueError, TypeError, OverflowError) as exc:
            _invalid(f"{name}.{column} cannot be restored as {dtype}: {exc}")
    return frame


def load_result(path):
    """Restore any result saved by ``result.save`` without fitting models.

    Args:
        path (str or pathlib.Path): Saved-result directory or its result.json file.
            JSON is authoritative; CSV files are convenience exports.

    Returns:
        EvaluationResult, AuditResult, SelectionResult, HierarchyEvaluationResult,
        or HierarchySelectionResult: Original result kind with tables and metadata.
        New saves preserve column order and dtypes, including empty tables.
        Earlier schema-1 files without column schemas infer dtypes from JSON.

    Raises:
        ValueError: Unsupported schema/kind or malformed JSON tables.
        OSError: The result file cannot be read.

    Examples:
        >>> import gauge_sc as gauge
        >>> evaluation = gauge.io.load_result("results/evaluation")
        >>> audit = gauge.audit(evaluation, reference="molecular", contrast="qc")
    """
    source = Path(path)
    if source.is_dir():
        source = source / "result.json"
    try:
        payload = json.loads(
            source.read_text(encoding="utf-8"), parse_constant=_nonfinite
        )
    except (json.JSONDecodeError, UnicodeError) as exc:
        _invalid(f"cannot decode JSON: {exc}")
    if not isinstance(payload, dict) or type(payload.get("schema_version")) is not int:
        _invalid("a versioned JSON object is required.")
    if payload["schema_version"] != 1:
        _invalid("unsupported schema_version; supported version is 1.")
    kind = payload.get("kind")
    if not isinstance(kind, str) or kind not in _RESULT_TYPES:
        _invalid("unknown result kind.")
    cls, names = _RESULT_TYPES[kind]
    metadata = payload.get("metadata")
    tables = payload.get("tables")
    schemas = payload.get("table_schema")
    if (
        not isinstance(metadata, dict)
        or not isinstance(tables, dict)
        or set(tables) != set(names)
    ):
        _invalid("metadata or required tables are missing.")
    if schemas is not None and (
        not isinstance(schemas, dict) or set(schemas) != set(names)
    ):
        _invalid("column schemas must cover every table.")
    frames = {
        name: _frame(tables[name], schemas[name] if schemas is not None else None, name)
        for name in names
    }
    for name, required in _REQUIRED[kind].items():
        if not required <= set(frames[name]):
            _invalid(f"table {name!r} is missing required columns.")
    if not _METADATA_KEYS[kind] <= set(metadata):
        _invalid("required result settings are missing from metadata.")
    if kind in ("evaluation", "hierarchy_evaluation") and not isinstance(
        metadata.get("views"), dict
    ):
        _invalid("evaluation metadata must include view definitions.")
    return cls(**frames, metadata=metadata)
