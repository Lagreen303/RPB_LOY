#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=60:00:00
#SBATCH --mem=150G
 
 
module load apptainer
 
output_dir="/iridisfs/ddnb/AIDA/LOY_velocyto"

cat $1 | while read a
do

apptainer run -H $PWD \
--bind /iridisfs/ddnb/AIDA/:/iridisfs/ddnb/AIDA/ \
--bind /iridisfs/home/lag1e24/images:/iridisfs/home/lag1e24/images/ \
--bind /iridisfs/ddnb/AIDA/LOY_cellranger:/iridisfs/ddnb/AIDA/LOY_cellranger \
--unsquash \
/iridisfs/home/lag1e24/images/datascience-notebook_latest.sif \
python3 run_call_LOY_velocyto.py \
$a $output_dir
done

 