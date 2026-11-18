"""Saved-result reuse and command-line contracts, including no-refit behavior."""

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier

import gauge_sc as ip
from gauge_sc.cli import main


@pytest.fixture
def inputs(tmp_path):
    codes = np.repeat(np.arange(4), 30)
    labels = np.array(["0", "001", "A", "B"])[codes]
    obs = pd.DataFrame(
        {
            "labels": labels,
            "coarse": np.where(codes < 2, "left", "right"),
            "depth": codes * 10.0 + 1.0,
        },
        index=[f"{i:04d}" for i in range(len(codes))],
    )
    data = ad.AnnData(codes[:, None].astype(float), obs=obs)
    data.obsm["X_pca"] = data.X.copy()
    source = tmp_path / "input.h5ad"
    data.write_h5ad(source)
    config = {
        "groupby": ["labels", "coarse"],
        "n_splits": 3,
        "views": {"molecular": {"obsm": "X_pca"}, "qc": {"obs": ["depth"]}},
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))
    return data, source, config, config_path


def _evaluation(inputs):
    data, _, config, _ = inputs
    return ip.evaluate(
        data, **config, classifier=DecisionTreeClassifier(random_state=0)
    )


def _assert_same(left, right):
    assert type(left) is type(right)
    assert left.metadata == right.metadata
    for name, value in vars(left).items():
        if name != "metadata":
            pd.testing.assert_frame_equal(value, getattr(right, name), check_exact=True)


def test_all_result_kinds_roundtrip_and_no_refit(inputs, tmp_path, monkeypatch):
    flat = _evaluation(inputs)
    data = inputs[0]
    hierarchy = ip.evaluate(
        data,
        hierarchy=["coarse", "labels"],
        views=inputs[2]["views"],
        n_splits=3,
        classifier=DecisionTreeClassifier(random_state=0),
    )
    audited = ip.audit(
        flat,
        groupby="labels",
        reference="molecular",
        contrast="qc",
        min_test_cells=100,
        maximum_mirror_ratio=1 / 3,
    )
    selected = ip.select(flat, view="molecular", delta=[0, 0.15])
    pruned = ip.select(hierarchy, view="molecular", delta=[0, 0.15])

    def no_fit(*args, **kwargs):
        raise AssertionError("Saved-result reuse must not fit a classifier")

    monkeypatch.setattr(DecisionTreeClassifier, "fit", no_fit)
    monkeypatch.setattr(RandomForestClassifier, "fit", no_fit)
    loaded = []
    for i, result in enumerate([flat, audited, selected, hierarchy, pruned]):
        directory = tmp_path / f"result-{i}"
        result.save(directory)
        restored = ip.io.load_result(directory)
        _assert_same(result, restored)
        _assert_same(result, ip.io.load_result(directory / "result.json"))
        loaded.append(restored)
        with pytest.raises(FileExistsError):
            result.save(directory)
    assert {"0", "001", "A", "B"} == set(
        loaded[0].groups.query("groupby == 'labels'").group
    )
    assert loaded[1].table.contrast.isna().all()
    assert not loaded[1].table.candidate.any()
    _assert_same(selected, ip.select(loaded[0], view="molecular", delta=[0, 0.15]))
    _assert_same(pruned, ip.select(loaded[3], view="molecular", delta=[0, 0.15]))
    _assert_same(
        audited,
        ip.audit(
            loaded[0],
            groupby="labels",
            reference="molecular",
            contrast="qc",
            min_test_cells=100,
            maximum_mirror_ratio=1 / 3,
        ),
    )
    # Restored writeback validates unchanged string IDs and group assignments.
    loaded[1].write_to_adata(data, key_added="restored_audit")
    loaded[4].write_to_adata(data, delta=0, key_added="restored_selection")


def test_empty_hierarchy_evidence_and_legacy_schema(inputs, tmp_path):
    hierarchy = ip.evaluate(
        inputs[0],
        hierarchy=["coarse", "labels"],
        views=inputs[2]["views"],
        min_child_cells=1000,
    )
    directory = tmp_path / "empty"
    hierarchy.save(directory)
    restored = ip.io.load_result(directory)
    _assert_same(hierarchy, restored)
    assert restored.predictions.empty
    _assert_same(
        ip.select(hierarchy, view="molecular", delta=0.1),
        ip.select(restored, view="molecular", delta=0.1),
    )
    flat = _evaluation(inputs)
    flat.save(tmp_path / "legacy")
    path = tmp_path / "legacy" / "result.json"
    payload = json.loads(path.read_text())
    del payload["table_schema"]
    path.write_text(json.dumps(payload))
    loaded = ip.io.load_result(path)
    _assert_same(
        ip.select(flat, view="molecular", delta=0.1),
        ip.select(loaded, view="molecular", delta=0.1),
    )


@pytest.mark.parametrize(
    "mutation",
    [
        "version",
        "kind",
        "table",
        "bool",
        "identity",
        "dtype",
        "columns",
        "required_column",
        "metadata",
        "nan",
    ],
)
def test_loader_rejects_malformed_results(inputs, tmp_path, mutation):
    result = _evaluation(inputs)
    result.save(tmp_path / "saved")
    source = tmp_path / "saved" / "result.json"
    payload = json.loads(source.read_text())
    if mutation == "version":
        payload["schema_version"] = 99
    elif mutation == "kind":
        payload["kind"] = "unknown"
    elif mutation == "table":
        del payload["tables"]["folds"]
    elif mutation == "bool":
        for record in payload["tables"]["table"]:
            record["made_up_flag"] = "false"
        payload["table_schema"]["table"]["made_up_flag"] = "bool"
    elif mutation == "identity":
        payload["tables"]["predictions"][0]["obs_name"] = 0
    elif mutation == "dtype":
        payload["table_schema"]["table"]["itca"] = "arbitrary"
    elif mutation == "columns":
        del payload["tables"]["groups"][0]["group"]
    elif mutation == "required_column":
        del payload["table_schema"]["table"]["itca"]
        for record in payload["tables"]["table"]:
            del record["itca"]
    elif mutation == "metadata":
        del payload["metadata"]["views"]
    elif mutation == "nan":
        payload["metadata"]["bad"] = float("nan")
    source.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="INVALID_RESULT_FILE"):
        ip.io.load_result(source)


def _run(args, capsys):
    code = main(args)
    captured = capsys.readouterr()
    if code == 0:
        assert not captured.err
        return json.loads(captured.out)
    assert not captured.out
    return code, json.loads(captured.err)


def test_cli_evaluate_audit_select_input_unchanged(
    inputs, tmp_path, capsys, monkeypatch
):
    _, source, _, config_path = inputs
    original_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    inventory = _run(["inspect", str(source)], capsys)
    assert inventory["n_cells"] == 120 and inventory["obsm"]["X_pca"] == [120, 1]
    directory = tmp_path / "evaluation"
    summary = _run(
        [
            "evaluate",
            str(source),
            "--config",
            str(config_path),
            "--output",
            str(directory),
        ],
        capsys,
    )
    assert summary["kind"] == "EvaluationResult"
    assert Path(summary["result_json"]).is_file()
    assert original_hash == hashlib.sha256(source.read_bytes()).hexdigest()
    evaluation = ip.io.load_result(directory)
    source.unlink()  # Audit/selection must work without the dataset.

    def no_fit(*args, **kwargs):
        raise AssertionError("CLI audit/select must not refit")

    monkeypatch.setattr(RandomForestClassifier, "fit", no_fit)
    audited_path = tmp_path / "audit"
    summary = _run(
        [
            "audit",
            str(directory),
            "--groupby",
            "labels",
            "--reference",
            "molecular",
            "--contrast",
            "qc",
            "--maximum-mirror-ratio",
            "0.3333333333333333",
            "--output",
            str(audited_path),
        ],
        capsys,
    )
    assert summary["screening_status"] == "no_threshold"
    _assert_same(
        ip.io.load_result(audited_path),
        ip.audit(
            evaluation,
            groupby="labels",
            reference="molecular",
            contrast="qc",
            maximum_mirror_ratio=1 / 3,
        ),
    )
    selected_path = tmp_path / "selection"
    summary = _run(
        [
            "select",
            str(directory),
            "--view",
            "molecular",
            "--delta",
            "0",
            "0.15",
            "--output",
            str(selected_path),
        ],
        capsys,
    )
    assert len(summary["path"]) == 2
    _assert_same(
        ip.io.load_result(selected_path),
        ip.select(evaluation, view="molecular", delta=[0, 0.15]),
    )
    # Output collision is checked before even trying to read the removed input.
    code, error = _run(
        [
            "evaluate",
            str(source),
            "--config",
            str(config_path),
            "--output",
            str(directory),
        ],
        capsys,
    )
    assert code == 2 and error["error"]["code"] == "OUTPUT_EXISTS"


def test_cli_hierarchy_selection(inputs, tmp_path, capsys):
    _, source, config, config_path = inputs
    config.pop("groupby")
    config["hierarchy"] = ["coarse", "labels"]
    config_path.write_text(json.dumps(config))
    directory = tmp_path / "tree"
    summary = _run(
        [
            "evaluate",
            str(source),
            "--config",
            str(config_path),
            "--output",
            str(directory),
        ],
        capsys,
    )
    assert summary["kind"] == "HierarchyEvaluationResult"
    selected = tmp_path / "pruned"
    summary = _run(
        [
            "select",
            str(directory),
            "--view",
            "molecular",
            "--population-ratio",
            "0.05",
            "0.1",
            "--output",
            str(selected),
        ],
        capsys,
    )
    assert summary["kind"] == "HierarchySelectionResult"
    _assert_same(
        ip.io.load_result(selected),
        ip.select(
            ip.io.load_result(directory), view="molecular", population_ratio=[0.05, 0.1]
        ),
    )


@pytest.mark.parametrize(
    "args,code",
    [
        (["audit"], "INVALID_ARGUMENT"),
        (["inspect", "missing.h5ad"], "IO_ERROR"),
        (
            [
                "select",
                "missing",
                "--view",
                "molecular",
                "--delta",
                "0.1",
                "--population-ratio",
                "0.1",
                "--output",
                "unused",
            ],
            "INVALID_ARGUMENT",
        ),
    ],
)
def test_cli_errors_are_machine_readable(args, code, capsys):
    status, error = _run(args, capsys)
    assert status == 2 and error["error"]["code"] == code


def test_cli_invalid_config_does_not_create_output(inputs, tmp_path, capsys):
    _, source, config, config_path = inputs
    config["classifier"] = "arbitrary-import"
    config_path.write_text(json.dumps(config))
    output = tmp_path / "bad"
    status, error = _run(
        [
            "evaluate",
            str(source),
            "--config",
            str(config_path),
            "--output",
            str(output),
        ],
        capsys,
    )
    assert status == 2 and error["error"]["code"] == "INVALID_CONFIG"
    assert not output.exists()


def test_installed_cli_entrypoints():
    # Installed console script and module entry point must agree on version.
    binary = Path(sys.executable).parent / ("gauge.exe" if os.name == "nt" else "gauge")
    for command in ([str(binary)], [sys.executable, "-m", "gauge_sc"]):
        result = subprocess.run(
            [*command, "--version"], capture_output=True, text=True, check=True
        )
        assert result.stdout.strip() == ip.__version__
        help_result = subprocess.run(
            [*command, "--help"], capture_output=True, text=True, check=True
        )
        assert all(
            name in help_result.stdout
            for name in ["inspect", "evaluate", "audit", "select"]
        )


def test_pandas3_string_schema_keeps_null_hierarchy_roots(inputs, tmp_path):
    result = ip.evaluate(
        inputs[0],
        hierarchy=["coarse", "labels"],
        views=inputs[2]["views"],
        n_splits=3,
        classifier=DecisionTreeClassifier(random_state=0),
    )
    result.save(tmp_path / "string_schema")
    source = tmp_path / "string_schema" / "result.json"
    payload = json.loads(source.read_text())
    payload["table_schema"]["nodes"]["parent_id"] = "str"
    source.write_text(json.dumps(payload))
    restored = ip.io.load_result(source)
    assert restored.nodes.parent_id.isna().sum() == 2
    _assert_same(
        ip.select(result, view="molecular", delta=0.1),
        ip.select(restored, view="molecular", delta=0.1),
    )
