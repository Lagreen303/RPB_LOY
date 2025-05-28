import pandas as pd
import numpy as np
import argparse
 
# Set up argument parser
parser = argparse.ArgumentParser(description="Merge LOY status from two files based on standardized barcodes.")
parser.add_argument('--file1', type=str, required=True, help="Path to the first input file (e.g., velocyto LOY status file).")
parser.add_argument('--file2', type=str, required=True, help="Path to the second input file (e.g., cellranger LOY status file).")
parser.add_argument('--output_file', type=str, required=True, help="Path to save the output merged file.")
 
# Parse the arguments
args = parser.parse_args()
 
# Load the input files
file1 = pd.read_csv(args.file1)
file2 = pd.read_csv(args.file2)
 
# Remove prefix and suffix in file1, and suffix in file2 to standardize barcodes
file1['barcode_standardized'] = file1['barcode'].str.replace(r'^.*:', '', regex=True).str.replace(r'x$', '', regex=True)
file2['barcode_standardized'] = file2['barcode'].str.replace(r'-\d+$', '', regex=True)
 
# Merge on the standardized barcode column
merged_data = pd.merge(file2, file1[['barcode_standardized', 'LOY_status']],
                       on='barcode_standardized', how='left', suffixes=('_file2', '_file1'))
 
# Fill NaN values with 0 for `LOY_status` columns to handle missing data
merged_data[['LOY_status_file2', 'LOY_status_file1']] = merged_data[['LOY_status_file2', 'LOY_status_file1']].fillna(0)
 
# Add `LOY_status_all` column: 1 if either LOY_status_file2 or LOY_status_file1 is 1
merged_data['LOY_status_all'] = np.where(
    (merged_data['LOY_status_file2'] == 1) & (merged_data['LOY_status_file1'] == 1),
    1, 0
)
 
# Save the merged result
merged_data.to_csv(args.output_file, index=False)
print(f"Output saved to {args.output_file}")
