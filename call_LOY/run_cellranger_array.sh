#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64 
#SBATCH --time=60:00:00
#SBATCH --partition=amd
#SBATCH --array=0-138

# Define the base directory containing all UUID directories
base_directory="/iridisfs/ddnb/AIDA"

# Define the path to the Cell Ranger executable
CELLRANGER_PATH="/iridisfs/home/lag1e24/images/cellranger-9.0.0/cellranger"

# Define the transcriptome reference
TRANSCRIPTOME_REF="/iridisfs/home/lag1e24/reference/refdata-gex-GRCh38-2024-A/"

# Read the list of directories from the file
mapfile -t dirs < ./directories_list

# Get the directory for this task
dir=${dirs[$SLURM_ARRAY_TASK_ID]}

# Extract the UUID from the directory name
uuid=$(basename "$dir")

# Define the FASTQ directory for the current UUID
FASTQ_DIR="$dir"

# Since the sample names are the same as the UUID, we can directly use the UUID
samples=${uuid}

# Create a unique output directory within the base directory
OUTDIR="${base_directory}/${uuid}/${uuid}_cellranger"

# Run Cell Ranger
$CELLRANGER_PATH count --id=${uuid}_cellranger \
 --sample=${samples} \
 --transcriptome=${TRANSCRIPTOME_REF} \
 --fastqs=${FASTQ_DIR} \
 --create-bam=true \
 --localcores=64 \
 --localmem=128 \

echo "Processed sample: $uuid"

