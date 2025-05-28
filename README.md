# LOY Project: Gene Set Enrichment and Classification Pipeline

This repository contains code and scripts used in a research project titled "Using Multi-omics to Characterise Clonal Expansion" focused on transcriptomic analysis and classification related to **loss of Y chromosome (LOY)** in single-cell RNA-seq data.

---

## Repository Structure

``
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

##  Overview
### 🔬 Biological Objective
To identify and classify transcriptional patterns associated with LOY using:
- Differential expression and gene set enrichment analysis (GSEA)
- Marker discovery
- Machine learning-based classification

---

## 🔧 Components

### `GSEA/heatmapcreation.r`
Scripts for gene set enrichment analysis and generating heatmaps. D

### `ML_models/`
Directory containing machine learning models and scripts used to classify LOY vs. non-LOY cells based on transcriptomic features.

### `call_LOY/`
Scripts for calling LOY using cellranger and velocyto.

### `find_markers.R`
Differential gene expression analysis R script to identify marker genes, using Seurat's `FindMarkers()` function. 

---

