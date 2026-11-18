# GAUGE frozen figure results

Redraw 41 quantitative panels from frozen results. The HTML overview shows
Figures 1–6. The individual PNG and PDF files use the paper's panel IDs. See
[Data sources and terms](DATA-SOURCES.md) for attribution and data access.
The archive includes the MIT license for its Python scripts.

## Run the replay

Download and unpack `gauge-frozen-results-v1.zip` from
[GAUGE Releases](https://github.com/zhanglabtools/gauge-sc/releases).
Run these commands in the unpacked directory with Python 3.12.

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python reproduction/fetch_cosmx_metadata.py \
  --output inputs/Run5642_S3_Quarter_metadata_file.csv
.venv/bin/python reproduction/render.py --bundle . --output rendered \
  --cosmx-metadata inputs/Run5642_S3_Quarter_metadata_file.csv
.venv/bin/python reproduction/validate.py --bundle . \
  --report rendered/scientific-checks.json
```

If the metadata fetch is unavailable, obtain `Quarter Brain.zip` from the
[Bruker download page](https://brukerspatialbiology.com/resources/coronal-hippocampus-and-cortex-basic-data-files/).
Extract `Run5642_S3_Quarter_metadata_file.csv` from it, place the CSV in
`inputs/`, then start with the render command. The renderer verifies the
metadata checksum. Replaying a selection such as `--figures 2 4` requires only
the bundled files. The `--figures` option uses paper figure numbers.

Open `rendered/index.html` for the figure overview. The `rendered/` directory
also contains the panel PNGs and PDFs and `scientific-checks.json`. The replay
starts from saved experiment outputs.
