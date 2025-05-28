library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(readr)
library(tibble)

# --- 0. Load gene annotation and exclude chrY genes ---
gene_info <- read_tsv("longest_transcripts_unique.txt")
y_genes <- gene_info %>%
  filter(chromosome == "chrY") %>%
  pull(gene)

# --- 1. Setup ---
input_dir <- "cell_type_DEGs/cell_type_DEGs_genename"
output_dir <- "GSEA_heatmaps"
dir.create(output_dir, showWarnings = FALSE)

DEG_files <- list.files(input_dir, pattern = "^differentially_expressed_genes_.*\\.csv$", full.names = TRUE)

# --- 2. Loop through each DEG file ---
for (file in DEG_files) {
  cell_type <- gsub("differentially_expressed_genes_|\\.csv", "", basename(file))
  
  # Load and filter DEGs
  deg_data <- read_csv(file, show_col_types = FALSE)
  deg_data <- deg_data %>%
    filter(!ensembl_id %in% y_genes, !is.na(avg_log2FC))
  
  # Prepare ranked gene list
  ranked_genes <- deg_data %>%
    mutate(gene = ensembl_id) %>%
    select(gene, avg_log2FC) %>%
    arrange(desc(avg_log2FC)) %>%
    mutate(entrez = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID",
                           keytype = "ENSEMBL", multiVals = "first")) %>%
    filter(!is.na(entrez))
  
  gene_list <- ranked_genes$avg_log2FC
  names(gene_list) <- ranked_genes$entrez

  # Run GSEA for BP only
  gsea_res <- tryCatch({
    gseGO(
      geneList = gene_list,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      keyType = "ENTREZID",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      verbose = FALSE
    )
  }, error = function(e) NULL)
  
  # Save top 10 results if GSEA returned anything
  if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
    gsea_df <- as.data.frame(gsea_res) %>%
      filter(p.adjust < 0.05) %>%
      arrange(desc(abs(NES))) %>%
      slice_head(n = 10)

    out_file <- file.path(output_dir, paste0(cell_type, "_GSEA_BP_top10.csv"))
    write_csv(gsea_df, out_file)
  }
}
