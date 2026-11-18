"""Small command-line interface over the public computational API."""

import argparse
import json
from pathlib import Path
import sys

import anndata as ad

from . import __version__, audit, evaluate, io, select


class _Parser(argparse.ArgumentParser):
    def error(self, message):
        raise ValueError(f"INVALID_ARGUMENT: {message}")


def _parser():
    parser = _Parser(
        prog="gauge",
        description="GAUGE: ITCA evaluation, view audit and ITCA+ selection.",
    )
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    inspect = commands.add_parser(
        "inspect", help="List AnnData inputs without fitting."
    )
    inspect.add_argument("input", type=Path, help="Input .h5ad file")
    ev = commands.add_parser("evaluate", help="Evaluate fixed groupings and views.")
    ev.add_argument("input", type=Path, help="Input .h5ad file")
    ev.add_argument(
        "--config", type=Path, required=True, help="JSON evaluate keyword arguments"
    )
    au = commands.add_parser(
        "audit", help="Audit a saved flat evaluation without refitting."
    )
    au.add_argument("result", type=Path, help="Evaluation directory or result.json")
    au.add_argument("--reference", required=True)
    au.add_argument("--contrast", required=True)
    au.add_argument("--groupby")
    au.add_argument("--min-test-cells", type=int, default=5)
    au.add_argument("--maximum-mirror-ratio", type=float)
    se = commands.add_parser(
        "select", help="Select from a saved flat/hierarchy evaluation."
    )
    se.add_argument("result", type=Path, help="Evaluation directory or result.json")
    se.add_argument("--view", required=True)
    penalty = se.add_mutually_exclusive_group(required=True)
    penalty.add_argument("--delta", type=float, nargs="+")
    penalty.add_argument("--population-ratio", type=float, nargs="+")
    for command in (ev, au, se):
        command.add_argument(
            "--output",
            type=Path,
            required=True,
            help="New result directory; never overwritten",
        )
    return parser


def _config(path):
    try:
        config = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, UnicodeError) as exc:
        raise ValueError(f"INVALID_CONFIG: cannot decode JSON: {exc}") from exc
    allowed = {
        "groupby",
        "hierarchy",
        "views",
        "n_splits",
        "random_state",
        "aggregation",
        "min_child_cells",
    }
    if not isinstance(config, dict) or set(config) - allowed or "views" not in config:
        raise ValueError(
            "INVALID_CONFIG: provide views and only supported evaluate keyword arguments."
        )
    return config


def _inspect(path):
    adata = ad.read_h5ad(path, backed="r")
    try:
        return {
            "kind": "input_inventory",
            "input": str(path.resolve()),
            "n_cells": adata.n_obs,
            "n_features": adata.n_vars,
            "obs_names_unique": bool(adata.obs_names.is_unique),
            "obs": {
                name: {
                    "dtype": str(series.dtype),
                    "n_unique": int(series.nunique()),
                    "n_missing": int(series.isna().sum()),
                }
                for name, series in adata.obs.items()
            },
            "obsm": {name: list(matrix.shape) for name, matrix in adata.obsm.items()},
            "layers": list(adata.layers),
        }
    finally:
        adata.file.close()


def _run(args):
    if args.command == "inspect":
        return _inspect(args.input)
    if args.output.exists():
        raise ValueError("OUTPUT_EXISTS: choose a new result directory.")
    if args.command == "evaluate":
        config = _config(args.config)
        result = evaluate(ad.read_h5ad(args.input), **config)
    else:
        evaluation = io.load_result(args.result)
        if args.command == "audit":
            result = audit(
                evaluation,
                reference=args.reference,
                contrast=args.contrast,
                groupby=args.groupby,
                min_test_cells=args.min_test_cells,
                maximum_mirror_ratio=args.maximum_mirror_ratio,
            )
        else:
            result = select(
                evaluation,
                view=args.view,
                delta=args.delta,
                population_ratio=args.population_ratio,
            )
    result.save(args.output)
    summary = {
        "kind": type(result).__name__,
        "output": str(args.output.resolve()),
        "result_json": str((args.output / "result.json").resolve()),
        "tables": {
            name: len(value)
            for name, value in vars(result).items()
            if name != "metadata"
        },
    }
    if args.command == "audit":
        summary.update(
            status=result.metadata["status"],
            screening_status=result.metadata["screening_status"],
            threshold=result.metadata["threshold"],
            screening=result.metadata["screening"],
        )
    if args.command == "select":
        summary["path"] = json.loads(result.path.to_json(orient="records"))
    return summary


def main(argv=None):
    """Return an exit status; successful commands emit one JSON object."""
    try:
        args = _parser().parse_args(argv)
        summary = _run(args)
    except (ValueError, TypeError, OSError) as exc:
        message = str(exc)
        prefix, separator, detail = message.partition(": ")
        code = (
            prefix
            if separator and prefix.isupper() and " " not in prefix
            else ("IO_ERROR" if isinstance(exc, OSError) else "INVALID_INPUT")
        )
        print(
            json.dumps(
                {
                    "error": {
                        "code": code,
                        "message": detail if code == prefix else message,
                    }
                }
            ),
            file=sys.stderr,
        )
        return 2
    print(json.dumps(summary, allow_nan=False))
    return 0
