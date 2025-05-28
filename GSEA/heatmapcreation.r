library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(circlize)
library(grid)

# --- 1. Setup ---
input_dir <- "GSEA_heatmaps"
output_dir <- "GSEA_heatmaps_by_lineage"
dir.create(output_dir, showWarnings = FALSE)

# Group GO terms into simplified categories
categorise_go_term <- function(description) {
  case_when(
    grepl("phagocytosis", description, ignore.case = TRUE) ~ "Phagocytosis",
    grepl("MHC|antigen|peptide", description, ignore.case = TRUE) ~ "Antigen Processing / MHC",
    grepl("complement activation", description, ignore.case = TRUE) ~ "Complement Activation",
    grepl("membrane invagination|vesicle", description, ignore.case = TRUE) ~ "Vesicle / Membrane Dynamics",
    grepl("biosynthetic|translation|ribosome", description, ignore.case = TRUE) ~ "Translation / Biosynthesis",
    grepl("interleukin|cytokine", description, ignore.case = TRUE) ~ "Cytokine Signaling",
    grepl("receptor|signaling|GPCR", description, ignore.case = TRUE) ~ "Receptor Signaling",
    grepl("T cell|B cell", description, ignore.case = TRUE) ~ "Lymphocyte Activation",
    grepl("ubiquitin|protein catabolic", description, ignore.case = TRUE) ~ "Protein Degradation",
    grepl("blood circulation|heart contraction", description, ignore.case = TRUE) ~ "Circulation / Heart",
    grepl("immune", description, ignore.case = TRUE) ~ "Immune Response",
    TRUE ~ description  # fallback to original description
  )
}


lineage_map <- c(
  "Treg" = "CD4_T", "CD8+_T_GZMB+" = "CD8_T", "CD8+_T_GZMK+" = "CD8_T",
  "CD4+_T_cyt" = "CD4_T", "gdT" = "Other_T", "dnT" = "Other_T", "CD4+_T" = "CD4_T",
  "CD4+_T_em" = "CD4_T", "CD8+_T" = "CD8_T", "CD4+_T_cm" = "CD4_T", 
  "CD4+_T_naive" = "CD4_T", "CD8+_T_naive" = "CD8_T", "MAIT" = "Other_T", "T" = "Other_T",
  "IGHMlo_memory_B" = "B", "IGHMhi_memory_B" = "B", "naive_B" = "B", 
  "atypical_B" = "B", "Plasma_B" = "B", "B" = "B",
  "CD16+_NK" = "NK", "NK" = "NK", "CD56+_NK" = "NK",
  "CD14+_Monocyte" = "Myeloid", "Monocyte" = "Myeloid", "CD16+_Monocyte" = "Myeloid",
  "cDC1" = "Myeloid", "cDC2" = "Myeloid", "cDC" = "Myeloid", "DC" = "Myeloid", "pDC" = "Myeloid"
)

lineage_colors <- c(
  "CD4_T" = "#1f77b4", "CD8_T" = "#2c7fb8", "Other_T" = "#6baed6",
  "B" = "#ff7f0e", "NK" = "#2ca02c", "Myeloid" = "#d62728", "Other" = "#aaaaaa"
)

top_n_per_lineage <- 50

# --- 2. Load and annotate all GSEA results ---
gsea_files <- list.files(input_dir, pattern = "_GSEA_BP_top10.csv$", full.names = TRUE)

gsea_combined <- lapply(gsea_files, function(file) {
  cell_type <- gsub("_GSEA_BP_top10\\.csv", "", basename(file))
  df <- read_csv(file, show_col_types = FALSE)
  df$CellType <- cell_type
  df$Lineage <- recode(cell_type, !!!lineage_map, .default = "Other")
  df
}) %>% bind_rows()

# --- 3. Build heatmap per lineage using top N unique GO terms ---
lineages <- unique(gsea_combined$Lineage)

for (lineage in lineages) {
  df_sub <- gsea_combined %>% filter(Lineage == lineage)
  if (nrow(df_sub) == 0) next

  # Select top N GO terms by max abs(NES)
  top_terms <- df_sub %>%
    group_by(Description) %>%
    summarise(max_abs_nes = max(abs(NES)), .groups = "drop") %>%
    arrange(desc(max_abs_nes)) %>%
    slice_head(n = top_n_per_lineage) %>%
    pull(Description)

  df_sub <- df_sub %>% filter(Description %in% top_terms)
df_sub <- df_sub %>%
  mutate(Description = categorise_go_term(Description)) %>%
  group_by(CellType, Description) %>%
  summarise(NES = mean(NES), .groups = "drop")



  # Construct NES matrix (GO term x CellType)
  mat <- df_sub %>%
    select(CellType, Description, NES) %>%
    pivot_wider(names_from = CellType, values_from = NES, values_fill = 0) %>%
    column_to_rownames("Description") %>%
    as.matrix()

  mat <- mat[, order(colnames(mat))]

  # Create and save heatmap
  color_scale <- colorRamp2(c(-3, 0, 3), c("#2166ac", "white", "#b2182b"))
  png_file <- file.path(output_dir, paste0("GSEA_BP_heatmap_", lineage, "_top", top_n_per_lineage, ".png"))

  png(png_file, width = 1200, height = 500)

Heatmap(
    mat,
    name = "NES",
    col = color_scale,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 20, fontface = "bold"),
    row_names_gp = gpar(fontsize = 16, fontface = "bold"),
    rect_gp = gpar(col = "white", lwd = 0.4),
    row_title = paste(lineage, "lineage"),
    row_title_gp = gpar(fontsize = 20, fontface = "bold"),
    row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 20)),
    column_names_max_height = unit(4, "cm")
  ) |> draw(padding = unit(c(10, 10, 40, 10), "mm")) 


  dev.off()
}
