# Prepare a tutorial-data release

Data versions are independent of Python package versions. The package catalog points to [tutorial-data-v1](https://github.com/zhanglabtools/gauge-sc/releases/tag/tutorial-data-v1). Keep each `.h5ad` as a separate attachment so readers download only what they need. Do not commit datasets or put them in a wheel.

## Prepared output

`scripts/package_tutorial_data.py` produces five validated h5ad attachments, `manifest.json`, `SHA256SUMS`, and an inventory/attribution README. It refuses existing output directories and checks every source asset against the preparation manifest before copying. Public repository configuration and upload are separate actions.

Publish the tutorial assets on [zhanglabtools/gauge-sc](https://github.com/zhanglabtools/gauge-sc). Prepare the release files and a candidate catalog in a new local staging directory. Use the existing [data sources and licenses](https://github.com/zhanglabtools/gauge-sc/blob/tutorial-data-v1/DATA-SOURCES.md) as the attribution record for v1.

```bash
python scripts/package_tutorial_data.py --input artifacts/tutorial-data-release-v1 --output artifacts/data-releases/tutorial-data-v1-configured --version tutorial-data-v1 --github-repo zhanglabtools/gauge-sc --catalog artifacts/tutorial-data-v1-catalog.json
```

After the assets are attached to the GitHub Release, verify that each public URL serves the pinned bytes before copying the staged catalog to `src/gauge_sc/tutorial_data.json`. See [the user guide](real-data-tutorials.md) for public downloads and local loading.

## Version and verification contract

- Review the five source/provenance records and applicable upstream redistribution terms before uploading.
- Publish the prepared bytes with their manifest and SHA256SUMS. Keep a published version immutable; changed bytes require a new data tag and package catalog.
- Install the newly built wheel in a clean environment and verify first download, cache reuse and offline loading for each dataset. Re-run the five real notebooks when their data or computation changes.
- Keep notebook outputs compact. Sphinx checks source fingerprints and renders stored real-data outputs. Core CI uses generated data.

The public loader verifies pinned hashes; a mutable online manifest is not trusted to replace them. No token appears in the package, notebook, catalog or documentation source.
