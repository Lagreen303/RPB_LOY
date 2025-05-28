import numpy as np
import pandas as pd 
import argparse
import velocyto as vcy #check this is installed in your container
import os
import re  # Import re module for regex

# Step 1: Define function to handle LOY calculation
def calculate_loy(loom_file, output_dir):
    # Step 2: Define the MSY region on chromosome Y
    msy_chromosome = 'Y'
    msy_start = 2781480
    msy_end = 56887902
 
    # Load the loom file using velocyto
    vlm = vcy.VelocytoLoom(loom_file)
 
    # Step 3: Filter MSY genes based on chromosome and start/end positions
    msy_gene_mask = (vlm.ra['Chromosome'] == msy_chromosome) & \
                    (vlm.ra['Start'] >= msy_start) & \
                    (vlm.ra['End'] <= msy_end)
 
    # Extract MSY gene expression (both spliced and unspliced)
    msy_spliced = vlm.S[msy_gene_mask, :]
    msy_unspliced = vlm.U[msy_gene_mask, :]
    msy_ampigous = vlm.A[msy_gene_mask, :]
 
    # Step 4: Combine spliced and unspliced read counts
    total_msy_expression = msy_spliced + msy_unspliced + msy_ampigous
 
    # Step 5: Identify cells with zero expression in MSY genes (LOY status)
    loy_status = (total_msy_expression.sum(axis=0) == 0).astype(int)  # 1 for LOY, 0 for no LOY
 
    # Step 6: Prepare barcodes and LOY status for CSV output
    barcodes = vlm.ca['CellID']
    loy_status_df = pd.DataFrame({
        'barcode': barcodes,
        'LOY_status': loy_status
    })
 
    # Step 7: Derive the output CSV file name based on the directory containing the loom file
    loom_dir = os.path.dirname(loom_file)
    sample_id = os.path.basename(loom_dir)
    output_csv = os.path.join(output_dir, f"{sample_id}_LOY_status.csv")  # Use sample_id in CSV name
 
    # Step 8: Save the LOY status to a CSV file
    loy_status_df.to_csv(output_csv, index=False)
    print(f"LOY status CSV saved to: {output_csv}")
 
# Step 9: Set up argument parsing for command-line usage
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate LOY status from a velocyto loom file")
    parser.add_argument("loom_file", type=str, help="Path to the input loom file")
    parser.add_argument("output_dir", type=str, help="Directory to save the output CSV file")
 
    args = parser.parse_args()
 
    # Call the function with the input loom file and output directory
    calculate_loy(args.loom_file, args.output_dir)