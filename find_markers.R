library(Seurat)
library(ggplot2)

# Load the Seurat object
seu <- readRDS("AIDA_male_LOY.rds")

# Unique cell types
cell_types <- unique(seu@meta.data$author_cell_type)

# Open a PDF device to save all heatmaps in one file
pdf("combined_heatmaps.pdf", width = 12, height = 8)

# Iterate over each cell type
for (cell_type in cell_types) {

  # Subset the Seurat object
  seu_subset <- subset(seu, subset = author_cell_type == as.character(cell_type))
  if (sum(seu_subset$LOY_status == 1) < 3) {next}

  # Renormalise data
  seu_subset <- NormalizeData(seu_subset)
  seu_subset <- FindVariableFeatures(seu_subset)
  seu_subset <- ScaleData(seu_subset)

  # Find differentially expressed genes
  markers <- FindMarkers(
    object = seu_subset,
    ident.1 = 1,
    ident.2 = 0,
    group.by = "LOY_status",
    test.use = "MAST"
  )

  # Save the results
output_file <- paste0("DGE_GSEA/differentially_expressed_genes_", cell_type, ".csv")
  write.csv(markers, output_file)
}