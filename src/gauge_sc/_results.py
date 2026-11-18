"""Serializable result containers; no fitted estimators are retained."""

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
import pandas as pd


def _records(frame):
    return frame.astype(object).where(pd.notna(frame), None).to_dict(orient="records")


@dataclass
class EvaluationResult:
    """Tables for candidates/views, groups, folds and aligned OOF predictions.

    Attributes:
        table (DataFrame): One row per candidate and view; scores use natural logs.
        groups (DataFrame): Group recall, information mass and contribution.
        folds (DataFrame): Group recall and support in each held-out fold.
        predictions (DataFrame): Observation identity, assignment, fold and prediction.
        metadata (dict): Versioned evaluation settings and view definitions.
    """

    table: pd.DataFrame
    groups: pd.DataFrame
    folds: pd.DataFrame
    predictions: pd.DataFrame
    metadata: dict

    def save(self, path):
        """Save JSON tables to a new directory (never overwrite)."""
        _save(self, path, "evaluation", ["table", "groups", "folds", "predictions"])


@dataclass
class SelectionResult:
    """Candidate-level ITCA+ scores and a path of selected grouping columns."""

    table: pd.DataFrame
    path: pd.DataFrame
    metadata: dict

    def save(self, path):
        """Save all candidates, selected path and provenance; never overwrite."""
        _save(self, path, "selection", ["table", "path"])


@dataclass
class HierarchyEvaluationResult:
    """Stored conditional scores, topology and per-parent evaluation evidence."""

    table: pd.DataFrame
    nodes: pd.DataFrame
    assignments: pd.DataFrame
    predictions: pd.DataFrame
    groups: pd.DataFrame
    folds: pd.DataFrame
    metadata: dict

    def save(self, path):
        _save(
            self,
            path,
            "hierarchy_evaluation",
            ["table", "nodes", "assignments", "predictions", "groups", "folds"],
        )


@dataclass
class HierarchySelectionResult:
    """Ancestor-closed split decisions and exhaustive terminal assignments."""

    table: pd.DataFrame
    path: pd.DataFrame
    nodes: pd.DataFrame
    splits: pd.DataFrame
    assignments: pd.DataFrame
    original_assignments: pd.DataFrame
    metadata: dict

    def save(self, path):
        _save(
            self,
            path,
            "hierarchy_selection",
            ["table", "path", "nodes", "splits", "assignments", "original_assignments"],
        )

    def write_to_adata(
        self,
        adata,
        *,
        population_ratio=None,
        delta=None,
        key_added="gauge_sc_selection",
    ):
        """Write one global pruning to the unchanged evaluated cells, explicitly."""
        from ._hierarchy import selected_delta

        chosen_delta = selected_delta(self, population_ratio, delta)
        if not isinstance(key_added, str) or not key_added:
            raise ValueError("INVALID_KEY: key_added must be nonempty.")
        if adata.is_view:
            raise ValueError("ANNDATA_VIEW: call .copy() before writeback.")
        original = self.original_assignments.set_index("obs_name").leaf_id
        if not adata.obs_names.is_unique or set(adata.obs_names) != set(original.index):
            raise ValueError("CELL_MISMATCH: writeback requires all evaluated cells.")
        levels = self.metadata["hierarchy"]
        if (
            any(level not in adata.obs for level in levels)
            or adata.obs[levels].isna().any().any()
        ):
            raise ValueError(
                "GROUPING_MISMATCH: original hierarchy is missing or changed."
            )
        paths = [
            json.dumps(list(row), ensure_ascii=False)
            for row in adata.obs[levels].astype(str).itertuples(index=False, name=None)
        ]
        if paths != original.reindex(adata.obs_names).tolist():
            raise ValueError("GROUPING_MISMATCH: original hierarchy changed.")
        if key_added in adata.obs or key_added in adata.uns:
            raise ValueError("OUTPUT_EXISTS: choose a different key_added.")
        selected = self.assignments.loc[
            self.assignments.delta == chosen_delta
        ].set_index("obs_name")
        adata.obs[key_added] = pd.Categorical(selected.label.reindex(adata.obs_names))
        adata.uns[key_added] = {
            "delta": chosen_delta,
            "metadata_json": json.dumps(self.metadata, allow_nan=False),
            "nodes_json": self.nodes.to_json(orient="records"),
            "assignments_json": selected.reset_index().to_json(orient="records"),
        }


@dataclass
class AuditResult:
    """Group-level contrasts and screening metadata, not significance tests."""

    table: pd.DataFrame
    assignments: pd.DataFrame
    metadata: dict

    def save(self, path):
        """Save tables and settings to a new directory (never overwrite)."""
        _save(self, path, "audit", ["table", "assignments"])

    def write_to_adata(self, adata, *, key_added="gauge_sc_audit"):
        """Explicitly add group scores to obs and a serializable snapshot to uns.

        Requires the same observation set and grouping assignments as evaluation.
        Row order may differ. Scores are group-level values broadcast to cells.
        Existing output keys are never overwritten. Returns None.
        """
        groupby = self.metadata["groupby"]
        if not isinstance(key_added, str) or not key_added:
            raise ValueError("INVALID_KEY: key_added must be a nonempty string.")
        if not adata.obs_names.is_unique or set(adata.obs_names) != set(
            self.assignments.obs_name
        ):
            raise ValueError(
                "CELL_MISMATCH: writeback requires the evaluated observation set."
            )
        if groupby not in adata.obs:
            raise ValueError(f"MISSING_GROUPING: {groupby!r} is absent from obs.")
        expected = self.assignments.set_index("obs_name").group.reindex(adata.obs_names)
        if adata.obs[groupby].isna().any() or not np.array_equal(
            adata.obs[groupby].astype(str).to_numpy(), expected.to_numpy()
        ):
            raise ValueError("GROUPING_MISMATCH: assignments changed after evaluation.")
        columns = ["contrast", "eligible"]
        if self.metadata["maximum_mirror_ratio"] is not None:
            columns.append("candidate")
        keys = [f"{key_added}_{col}" for col in columns]
        if key_added in adata.uns or any(key in adata.obs for key in keys):
            raise ValueError("OUTPUT_EXISTS: choose a different key_added.")
        if adata.is_view:
            raise ValueError("ANNDATA_VIEW: call .copy() before explicit writeback.")
        lookup = self.table.set_index("group")
        for col, key in zip(columns, keys):
            adata.obs[key] = expected.map(lookup[col]).to_numpy()
        adata.uns[key_added] = {
            "metadata_json": json.dumps(self.metadata, allow_nan=False),
            "table_json": self.table.to_json(orient="records"),
            "score_level": "group; broadcast to observations",
        }


def _save(result, path, kind, tables):
    destination = Path(path)
    payload = {
        "schema_version": 1,
        "kind": kind,
        "metadata": result.metadata,
        "tables": {name: _records(getattr(result, name)) for name in tables},
        "table_schema": {
            name: {
                column: str(dtype)
                for column, dtype in getattr(result, name).dtypes.items()
            }
            for name in tables
        },
    }
    encoded = json.dumps(payload, indent=2, allow_nan=False)
    destination.mkdir(parents=True, exist_ok=False)
    (destination / "result.json").write_text(encoded + "\n", encoding="utf-8")
    for name in tables:
        getattr(result, name).to_csv(destination / f"{name}.csv", index=False)
