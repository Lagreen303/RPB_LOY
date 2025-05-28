#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64 
#SBATCH --time=60:00:00
#SBATCH --mem=150G
#SBATCH --partition=amd
#SBATCH --array=0-138

# Read the list of directories from the file
mapfile -t dirs < ./directories_list

# Get the directory for this task
sample=${dirs[$SLURM_ARRAY_TASK_ID]}

echo $sample

module load apptainer

# Bind my applications and the shared data spaces (or wherever the data is stored)
export run_velocyto="apptainer exec -H $PWD \
--bind /iridisfs/ddnb/AIDA/:/iridisfs/ddnb/AIDA/ \
--bind /iridisfs/home/lag1e24/images/ref_mask/:/iridisfs/home/lag1e24/images/ref_mask/ \
--bind /iridisfs/home/lag1e24/reference/refdata-gex-GRCh38-2024-A/genes/ \
--bind /home/lag1e24/ddnb/data/AIDA/analysis/:/home/lag1e24/ddnb/data/AIDA/analysis/ \
--unsquash \
/iridisfs/home/lag1e24/images/velocyto_0.17.sif \
velocyto"

filtered_barcodes=/iridisfs/ddnb/AIDA/${sample}_cellranger/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
outdir=/iridisfs/ddnb/AIDA/${sample}_velocyto
bamfile="/iridisfs/ddnb/AIDA/${sample}_cellranger/outs/possorted_genome_bam.bam" #check this path and change if needed (change analysis to AIDA/sample/sample_cellranger)
repeat_msk="/iridisfs/home/lag1e24/images/ref_mask/repeat_msk.gtf"
GTFFILE=/iridisfs/ddnb/AIDA/genes.gtf #check this path and change if needed

echo GTFFILE: $GTFFILE

$run_velocyto run -o $outdir -b $filtered_barcodes -m $repeat_msk -@ $SLURM_CPUS_PER_TASK $bamfile $GTFFILE