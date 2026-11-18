"""Prepare GitHub Release attachments without committing or uploading datasets."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--version", default="tutorial-data-v1")
    parser.add_argument(
        "--github-repo",
        help="Confirmed owner/repo; omit while destination is undecided.",
    )
    parser.add_argument(
        "--catalog",
        type=Path,
        help="Explicitly update the small package catalog at this path.",
    )
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", args.version):
        parser.error("version must be a simple release tag")
    if args.github_repo and not re.fullmatch(
        r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", args.github_repo
    ):
        parser.error("github-repo must be owner/repo")
    catalog = json.loads((args.input / "manifest.json").read_text())
    catalog["version"] = args.version
    catalog["schema_version"] = 1
    catalog["release_destination"] = args.github_repo
    catalog.pop("publication_status", None)
    names = []
    for record in catalog["datasets"].values():
        name = record["filename"]
        if Path(name).name != name or not name.endswith(".h5ad") or name in names:
            raise ValueError("Invalid or duplicated asset filename")
        names.append(name)
        path = args.input / name
        if path.stat().st_size != record["bytes"] or digest(path) != record["sha256"]:
            raise ValueError(f"Asset does not match preparation manifest: {name}")
        record["url"] = (
            f"https://github.com/{args.github_repo}/releases/download/{args.version}/{name}"
            if args.github_repo
            else None
        )
    args.output.mkdir(parents=True, exist_ok=False)
    for name in names:
        shutil.copy2(args.input / name, args.output / name)
    text = json.dumps(catalog, indent=2) + "\n"
    (args.output / "manifest.json").write_text(text)
    (args.output / "SHA256SUMS").write_text(
        "".join(
            f"{digest(args.output / name)}  {name}\n"
            for name in [*names, "manifest.json"]
        )
    )
    rows = [
        f"| {key} | {row['n_obs']:,} | {row['n_vars']:,} | {row['bytes']/1e6:.2f} |"
        for key, row in catalog["datasets"].items()
    ]
    (args.output / "README.md").write_text(
        f"# ITCA+ tutorial data: {args.version}\n\n"
        "Versioned release attachments. The preparation script itself never uploads files.\n\n"
        "| Dataset | Cells | Features | MB |\n|---|---:|---:|---:|\n"
        + "\n".join(rows)
        + "\n\n"
        "Distribute the h5ad files as separate Release attachments, with manifest.json and SHA256SUMS. "
        "Do not git-add them or include them in Python wheels. Each tutorial loads only its own file. "
        "Freeze these checksums after publication; changed bytes need a new data version.\n\n"
        "The manifest records source identities, source checksums, sampling policies and expression scales. "
        "Allen Cortex contains normalized log1p expression over selected HVGs; ABC HPF contains supplied log2 expression. "
        "Do not log-transform these a second time. Ding mixture withheld species annotations are post-audit readouts.\n\n"
        "Keep upstream attribution and applicable redistribution terms with the public data release. "
        "This preparation script does not assign a new license to upstream datasets.\n"
    )
    if args.catalog:
        args.catalog.write_text(text)
    print(
        f"Prepared {len(names)} assets ({sum(r['bytes'] for r in catalog['datasets'].values())/1e6:.2f} MB): {args.output}"
    )
    if not args.github_repo:
        print("Public download URLs are unset until an owner/repo is supplied.")


if __name__ == "__main__":
    main()
