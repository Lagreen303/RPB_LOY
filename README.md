# LOY Project: Gene Set Enrichment and Classification Pipeline

This repository contains code and scripts used in a research project focused on transcriptomic analysis and classification related to **loss of Y chromosome (LOY)** in single-cell RNA-seq data.

---

## 📁 Repository Structure

```
.
├── GSEA/
│   └── heatmapcreation.r
├── ML_models/
│   └── [ML model files here]
├── call_LOY/
│   └── merge_cellranger_velocyto_LOY_calls.py
├── find_markers.R
└── README.md
```

---

## 📊 Overview

### 🔬 Biological Objective
To identify and classify transcriptional patterns associated with LOY using:
- Differential expression and gene set enrichment analysis (GSEA)
- Marker discovery
- Machine learning-based classification

---

## 🔧 Components

### `GSEA/heatmapcreation.r`
Script for generating heatmaps from gene set enrichment results. Designed to visualize lineage-specific or group-specific GO:BP term enrichments.

### `ML_models/`
Directory containing machine learning models and scripts used to classify LOY vs. non-LOY cells based on transcriptomic features.

### `call_LOY/merge_cellranger_velocyto_LOY_calls.py`
Python script for merging Cell Ranger and velocyto outputs with LOY calls (likely derived from genomic or computational analysis), preparing integrated metadata for downstream analysis.

### `find_markers.R`
R script to identify marker genes, possibly using Seurat's `FindMarkers()` function. Useful for pinpointing key transcriptional differences between cell groups.

---

## 📦 Dependencies

### R
- Seurat
- clusterProfiler
- ggplot2
- ComplexHeatmap

### Python
- pandas
- scanpy or anndata
- numpy

---

## ▶️ Usage

You may need to run the following in order:
1. `call_LOY/merge_cellranger_velocyto_LOY_calls.py` to prepare metadata
2. `find_markers.R` to identify DEGs
3. GSEA pipeline including `GSEA/heatmapcreation.r`
4. Run classification models in `ML_models/`

---

## 📌 Notes

- Ensure all input files (e.g., Seurat/AnnData objects) are formatted consistently.
- Directory names are lowercase for clarity, but can be modified depending on internal standards.
- Consider setting environment variables or config files for paths and parameters.

---

## 🧪 Citation / Acknowledgement

This code is part of an ongoing research project. Please contact the maintainer if you intend to reuse or adapt the workflows.
