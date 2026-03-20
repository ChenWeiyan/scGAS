# scGAS <img src="man/figures/logo.png" align="right" height="139" />

**Single-Cell Gene Activation Potential Inference via a Gene Regulatory Reference Map**

[![R-CMD-check](https://github.com/yourname/scGAS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yourname/scGAS/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

`scGAS` infers **Gene Activation Potential (GAS)** from single-cell ATAC-seq
data by combining:

- A **reference map** of cCRE-Gene associations built from 167 paired ENCODE
  bulk DNase-seq / RNA-seq samples.
- **Per-gene Lasso regression models** that quantify the contribution of each
  cis-regulatory element to gene activation.
- A **Metacell + network propagation** strategy that overcomes the extreme
  sparsity of scATAC-seq data.
- A **Chromatin Potential Field** that predicts differentiation trajectories by
  linking chromatin state to RNA expression.

The method is described in:

> 

---

## Installation

```r
# Install development version from GitHub
remotes::install_github("ChenWeiyan/scGAS")

# Required Bioconductor packages
BiocManager::install(c("Signac", "EnsDb.Hsapiens.v75", "SCAVENGE"))
```

---

## Quick Start

```r
library(scGAS)
library(EnsDb.Hsapiens.v75)

# 1. Gene annotations
annotation <- Signac::GetGRangesFromEnsDb(EnsDb.Hsapiens.v75)
seqlevelsStyle(annotation) <- "UCSC"

# 2. Preprocess scATAC-seq
obj <- scgas_preprocess(
  fragment_paths     = c("rep1_fragments.tsv.gz", "rep2_fragments.tsv.gz"),
  cre_bed_path       = "hg19_500bp_CRE.bed",
  annotation         = annotation,
  cluster_resolution = 1.4,
  n_cores            = 8
)

# 3. Train per-gene Lasso models
models <- scgas_train_models(
  seurat_obj           = obj,
  sig_associations     = readRDS("ENCODE_167Tissue_SigAss_pVal.05.rds"),
  encode_rna_lognorm   = readRDS("ENCODE_167Tissue_RNA_LogNorm.rds"),
  encode_dnase_lognorm = readRDS("ENCODE_167Tissue_CRE_LogNorm.rds"),
  cre_bed_path         = "hg19_500bp_CRE.bed",
  n_cores              = 8
)

# 4. Compute single-cell GAS
scgas_mat <- scgas_compute(
  seurat_obj     = obj,
  trained_models = models,
  n_cores        = 8
)

# 5. Add to Seurat object and visualise
obj <- scgas_add_assay(obj, scgas_mat)

# 6. Chromatin Potential Field
cpf <- scgas_chromatin_potential(obj, point_colour = "celltype")
print(cpf$plot)
```

---

## Reference data

The following pre-computed files are required and can be downloaded from
[Zenodo](https://zenodo.org) *(10.5281/zenodo.19134899)*:

| File | Description |
|------|-------------|
| `hg19_500bp_CRE.bed` | 1.8 M reference cCRE bins (hg19) |
| `ENCODE_167Tissue_SigAss_pVal.05.rds` | Significant cCRE-Gene pairs |
| `ENCODE_167Tissue_RNA_LogNorm.rds` | ENCODE bulk RNA-seq matrix |
| `ENCODE_167Tissue_CRE_LogNorm.rds` | ENCODE bulk DNase-seq matrix |

---

## Main functions

| Function | Description |
|----------|-------------|
| `scgas_preprocess()` | Fragment → cCRE count matrix, QC, LSI, high-res clustering |
| `scgas_train_models()` | Metacell aggregation + per-gene Lasso training |
| `scgas_compute()` | Single-cell GAS via initialisation + network propagation |
| `scgas_add_assay()` | Embed scGAS matrix into Seurat object |
| `scgas_chromatin_potential()` | Build and plot the Chromatin Potential Field |

---

## Citation

If you use scGAS please cite:

```
```

---

## License

MIT © 2026 Weiyan Chen
