# Tutorial data sources

These five AnnData files are prepared for the GAUGE tutorials. The upstream license applies to each file individually.

| File | Cells × features | Upstream terms |
| --- | ---: | --- |
| `abc_hpf.h5ad` | 3,000 × 500 | [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/) |
| `allen_cortex.h5ad` | 3,000 × 2,000 | [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/) |
| `vasculature.h5ad` | 2,500 × 26,666 | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) |
| `ding_pbmc.h5ad` | 2,000 × 33,694 | [Single Cell Portal public-study terms](https://singlecell.broadinstitute.org/single_cell/terms_of_service) |
| `ding_mixture.h5ad` | 3,159 × 15,796 | [Single Cell Portal public-study terms](https://singlecell.broadinstitute.org/single_cell/terms_of_service) |

## Allen Institute whole mouse brain atlas

Both Allen files are derived from the Allen Institute for Brain Science's [Allen Brain Cell Atlas](https://alleninstitute.github.io/abc_atlas_access/). Cite Yao et al., [A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain](https://doi.org/10.1038/s41586-023-06812-z).

- **`abc_hpf.h5ad`** uses the [MERFISH expression data](https://alleninstitute.github.io/abc_atlas_access/descriptions/MERFISH-C57BL6J-638850.html) and [CCF-mapped coordinates](https://alleninstitute.github.io/abc_atlas_access/descriptions/MERFISH-C57BL6J-638850-CCF.html) for specimen C57BL6J-638850. The tutorial uses supplied log2 expression and retains a seeded sample from DG, CA1, CA2 and CA3 in the 6.8–7.0 mm CCF slice.
- **`allen_cortex.h5ad`** uses the [WMB-10X transcriptomes](https://alleninstitute.github.io/abc_atlas_access/descriptions/WMB-10X.html) and [WMB taxonomy](https://alleninstitute.github.io/abc_atlas_access/descriptions/WMB-taxonomy.html). The tutorial samples 3,000 cells from a 20,000-cell isocortex reference drawn from WMB-10Xv2 and WMB-10Xv3 and retains 2,000 highly variable genes after library-size normalization and log1p transformation.

## Wälchli et al. human brain vasculature atlas

**`vasculature.h5ad`** is derived from [Wälchli and Ghobrial's Figshare archive of FACS-sorted brain endothelial cells](https://doi.org/10.6084/m9.figshare.25151738), file ID 44431100, accompanying [Wälchli et al. (2024)](https://doi.org/10.1038/s41586-024-07493-y). The Seurat archive was sampled to 150,000 cells by pathology, donor and endothelial cluster, then exported as an AnnData count matrix. A seeded sample of 50,000 cells was masked for workflow comparison. The tutorial samples 2,500 cells from that input, retaining the counts and two count-based quality-control fields.

## Ding et al. single-cell sequencing comparison

Both inputs derive from [Ding et al. (2020)](https://doi.org/10.1038/s41587-020-0465-8), [GEO GSE132044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132044). The authors' corresponding Single Cell Portal studies are [SCP424 (PBMC)](https://singlecell.broadinstitute.org/single_cell/study/SCP424/single-cell-comparison-pbmc-data) and [SCP426 (mixture)](https://singlecell.broadinstitute.org/single_cell/study/SCP426/single-cell-comparison-mixture-data). The portal terms provide for public redistribution and reuse of data that contributors submit as public studies. These subsets are shared for academic research and teaching.

- **`ding_pbmc.h5ad`** includes 1,000 cells each from the 10x Chromium v3 and Drop-seq arms of PBMC1, sampled from cells with an assigned type without balancing types.
- **`ding_mixture.h5ad`** includes all 3,159 cells in the Mixture1 10x Chromium v2 sample. Counts are summed over matched human and mouse gene symbols. The mixed-species label is retained for post-audit validation.
