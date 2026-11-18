# Reproduce Figure 4b–f from source data

Figure 4 compares species-blind audits on two Ding et al. human–mouse mixtures (GSE132044). This workflow verifies a pinned source asset, redraws the five quantitative panels and checks the paper's displayed summaries.

## Run the standalone replay

Download and unpack `gauge-frozen-results-v1.zip` from [GAUGE Releases](https://github.com/zhanglabtools/gauge-sc/releases). Run these commands from the unpacked directory. Figures 4b–f use the bundled tables.

```bash
python -m pip install -r requirements.txt
python reproduction/render.py --bundle . --figures 4 --output rendered-figure4
python reproduction/validate.py --bundle . --figures 4 \
  --report rendered-figure4/scientific-checks.json
```

Open `rendered-figure4/index.html` for the preview. The same directory contains the five PNG/PDF pairs and `scientific-checks.json`.

| Panel | Result |
| --- | --- |
| 4b | Cluster contrast scores in both mixtures |
| 4c | QC and mixed-species fractions in Mixture 1 |
| 4d | QC and mixed-species fractions in Mixture 2 |
| 4e | ROC curves for multiplet-enriched clusters |
| 4f | AUROCs under matched randomization |

For a one-sample application of the public API, use the [real-data mixture audit tutorial](../notebooks/real_data/11_mixture_audit.ipynb).
