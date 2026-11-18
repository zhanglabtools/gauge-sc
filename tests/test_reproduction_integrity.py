"""Reject incomplete or corrupted paper source assets before plotting."""

import importlib
import json
from pathlib import Path

import pytest


@pytest.fixture
def asset(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).parents[1] / "reproduction"))
    render = importlib.import_module("render")
    bundle = tmp_path / "asset"
    (bundle / "reference").mkdir(parents=True)
    (bundle / "source_data").mkdir()
    authority = bundle / "reference/GAUGE-Main-Figures.pptx"
    authority.write_bytes(b"frozen-authority")
    source = bundle / "source_data/data.csv"
    source.write_text("value\n1\n")
    panels = []
    for figure, count in render.PANEL_COUNTS.items():
        for index in range(count):
            panel = f"fig{figure}{chr(97 + index)}"
            tables = [] if panel in render.SCHEMATICS else ["data.csv"]
            if figure == 4 and tables:
                panel_source = bundle / "source_data" / f"{panel}.csv"
                panel_source.write_bytes(source.read_bytes())
                tables = [panel_source.name]
            panels.append(
                dict(
                    id=panel,
                    figure=figure,
                    panel=chr(97 + index),
                    tables=tables,
                )
            )
    manifest = dict(
        schema_version=1,
        panels=panels,
        authority={"filename": authority.name, "sha256": render.digest(authority)},
        files={
            p.relative_to(bundle).as_posix(): {
                "bytes": p.stat().st_size,
                "sha256": render.digest(p),
            }
            for p in [authority, *(bundle / "source_data").glob("*.csv")]
        },
    )
    (bundle / "manifest.json").write_text(json.dumps(manifest))
    return bundle, manifest, render


def test_changed_source_is_rejected(asset):
    bundle, _, render = asset
    render.verify(bundle)
    (bundle / "source_data/data.csv").write_text("value\n2\n")
    with pytest.raises(ValueError, match="Checksum mismatch"):
        render.verify(bundle)


def test_missing_panel_is_rejected(asset):
    bundle, manifest, render = asset
    manifest["panels"].pop()
    (bundle / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="Panel inventory"):
        render.verify(bundle)


def test_unpinned_table_is_rejected(asset):
    bundle, manifest, render = asset
    manifest["files"].pop("source_data/data.csv")
    (bundle / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="not pinned"):
        render.verify(bundle)


def test_figure4_package_verifies_without_other_figures(asset, tmp_path):
    bundle, manifest, render = asset
    packager = importlib.import_module("package_figure4")
    manifest["limitations"] = [
        "Figure 2h: agreement pending.",
        "Figure 4e: check the upstream environment.",
        "Shared redistribution review pending.",
    ]
    (bundle / "manifest.json").write_text(json.dumps(manifest))
    output = tmp_path / "figure4-asset"
    packager.package(bundle, output)
    scoped = render.verify(output)
    assert scoped["figures"] == [4]
    assert [p["id"] for p in scoped["panels"]] == [f"fig4{p}" for p in "bcdef"]
    assert scoped["authority"] == manifest["authority"]
    assert not list(output.rglob("*.pptx"))
    assert sorted(p.name for p in (output / "source_data").glob("*.csv")) == [
        f"fig4{p}.csv" for p in "bcdef"
    ]
    assert scoped["limitations"] == manifest["limitations"][1:]
    assert scoped["source_asset_manifest_sha256"] == render.digest(
        bundle / "manifest.json"
    )

    scoped["panels"].pop()
    (output / "manifest.json").write_text(json.dumps(scoped))
    with pytest.raises(ValueError, match="Panel inventory"):
        render.verify(output)
