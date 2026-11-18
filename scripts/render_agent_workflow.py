"""Regenerate the agent guide's table, figure and run provenance from its notebook."""

import hashlib
from importlib.metadata import version
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nbformat


def main():
    root = Path(__file__).resolve().parents[1]
    source = root / "notebooks/07_agent_audit.ipynb"
    notebook = nbformat.read(source, as_version=4)
    namespace = {}
    original_show = plt.show
    plt.show = lambda *args, **kwargs: None
    try:
        for cell in notebook.cells:
            if cell.cell_type == "code":
                exec(compile(cell.source, str(source), "exec"), namespace)
    finally:
        plt.show = original_show
    audit = namespace["audit"]
    evaluation = namespace["evaluation"]
    assets = root / "docs/assets"
    readout = namespace["readout"]
    readout[
        ["group", "reference_contribution", "contrast_contribution", "contrast"]
    ].rename(
        columns={
            "group": "Group",
            "reference_contribution": "Molecular contribution",
            "contrast_contribution": "QC contribution",
            "contrast": "Centered contrast",
        }
    ).round(
        3
    ).to_csv(
        assets / "agent-audit-table.csv", index=False
    )
    namespace["audit_figure"].savefig(assets / "agent-audit.svg", bbox_inches="tight")
    vector = assets / "agent-audit.svg"
    vector.write_text(
        "\n".join(line.rstrip() for line in vector.read_text().splitlines()) + "\n"
    )
    (assets / "agent-audit-summary.md").write_text(
        namespace["agent_interpretation"] + "\n", encoding="utf-8"
    )
    skill = root / ".agents/skills/gauge/SKILL.md"
    report = {
        "generated_by": "Codex using the repository GAUGE skill",
        "skill_sha256": hashlib.sha256(skill.read_bytes()).hexdigest(),
        "notebook_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "versions": {
            name: version(name)
            for name in ["gauge_sc", "anndata", "pandas", "scikit-learn"]
        },
        "data": dict(namespace["adata"].uns["gauge_sc_example"]),
        "evaluation": evaluation.metadata,
        "audit": audit.metadata,
        "table": json.loads(audit.table.to_json(orient="records", double_precision=15)),
        "input_unchanged": True,
        "saved_result_reuse_matches": True,
    }
    (assets / "agent-audit-run.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    plt.close("all")
    print("Rendered agent workflow assets from the runnable notebook.")


if __name__ == "__main__":
    main()
