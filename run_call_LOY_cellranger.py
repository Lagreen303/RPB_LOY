import gzip
import pandas as pd
import scanpy as sc
import argparse
import os
import re

# Step 1: Setup argument parser for command-line inputs
def parse_args():
    parser = argparse.ArgumentParser(description="Detect LOY genes from Cell Ranger data.")
    parser.add_argument("cellranger_path", type=str, help="Path to the filtered_feature_bc_matrix directory from Cell Ranger output")
    parser.add_argument("gtf_file", type=str, help="Path to the GTF file (genes.gtf or genes.gtf.gz)")
    parser.add_argument("output_dir", type=str, help="Directory to save the output CSV file")
    return parser.parse_args()

# Step 2: Define function to detect LOY genes
def detect_loy(cellranger_path, gtf_file, output_dir):
    # Define the MSY region on chromosome Y
    msy_chromosome = 'chrY'
    msy_start = 2781480
    msy_end = 56887902

    # Step 3: Determine whether the GTF file is gzipped or not, and open it accordingly
    if gtf_file.endswith('.gz'):
        open_func = gzip.open
        mode = "rt"
    else:
        open_func = open
        mode = "r"

    # Load the GTF file and filter for genes in the MSY region
    msy_genes = []

    with open_func(gtf_file, mode) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] == "gene" and fields[0] == msy_chromosome:
                gene_start = int(fields[3])
                gene_end = int(fields[4])
                if gene_start >= msy_start and gene_end <= msy_end:
                    # Extract gene_id from attributes column (last column in the GTF)
                    gene_info = fields[-1]
                    gene_name = gene_info.split('gene_id "')[1].split('"')[0]
                    msy_genes.append(gene_name)

    # Convert MSY genes list to a set for fast lookup
    msy_gene_set = set(msy_genes)
    print(f"Total MSY genes identified: {len(msy_gene_set)}")
    print(msy_gene_set)
    # Step 4: Load Cell Ranger output using scanpy
    adata = sc.read_10x_mtx(cellranger_path, var_names='gene_ids', cache=True)
    print(f"Data shape: {adata.shape}")  # (n_cells, n_genes)

    # Step 5: Filter for MSY genes
    msy_gene_mask = adata.var_names.isin(msy_gene_set)

    # Subset the AnnData object to only include MSY genes
    adata_msy = adata[:, msy_gene_mask].copy()
    print(f"Shape of MSY data: {adata_msy.shape}")  # Should be (n_cells, n_MSY_genes)

    # Step 6: Find cells with zero expression in the MSY region
    no_msy_expression_cells = (adata_msy.X.sum(axis=1) == 0).A1  # Boolean array where True = no MSY expression

    # Get barcodes (cell names) for cells with no MSY expression
    no_msy_expression_barcodes = adata.obs_names[no_msy_expression_cells]
    print(f"Number of cells with no expression in MSY region: {len(no_msy_expression_barcodes)}")

    # Step 7: Assign LOY Status and Generate CSV File
    adata.obs['LOY_status'] = 0  # Default: No LOY
    adata.obs.loc[no_msy_expression_barcodes, 'LOY_status'] = 1  # Mark LOY cells as 1

    # Generate a DataFrame with barcode and LOY status
    loy_status_df = pd.DataFrame({
        'barcode': adata.obs_names,
        'LOY_status': adata.obs['LOY_status'].values
    })

    # Step 8: Generate the output CSV path using the specified output directory and sample ID
    match = re.search(r'/iridisfs/ddnb/AIDA/([^/]+)/outs/filtered_feature_bc_matrix', cellranger_path)
    if match:
        sample_id = match.group(1)
    else:
        sample_id = os.path.basename(os.path.normpath(cellranger_path))
    output_csv_path = os.path.join(output_dir, f"{sample_id}_LOY_status.csv")

    # Save the DataFrame as a CSV file
    loy_status_df.to_csv(output_csv_path, index=False)
    print(f"LOY status CSV saved to: {output_csv_path}")

# Step 9: Main function to execute the script
if __name__ == "__main__":
    args = parse_args()
    detect_loy(args.cellranger_path, args.gtf_file, args.output_dir)