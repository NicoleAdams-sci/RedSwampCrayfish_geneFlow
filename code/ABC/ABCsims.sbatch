#!/bin/bash
#SBATCH --time=23:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --job-name=rscABC
##SBATCH --error=../QSTAT/ABC.%J.stderr
#SBATCH --output=../QSTAT/ABC.%J.stdout
#SBATCH --constraint="[intel18|amd20]"

# Usage:  
# run from /mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SHELL

module purge
module load GCC/12.2.0 OpenMPI/4.1.4-GCC-12.2.0 R/4.2.2-foss-2022b ICU/72.1-GCCcore-12.2.0

DIR="/mnt/ufs18/rs-012/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/"

source $DIR/setup_R.sh

#Rscript ../RSC_ABCsims_10-2-2023.R ${SLURM_ARRAY_TASK_ID} 100 ${USER}
$DIR/SRC/R-4.1.2-install/bin/Rscript $DIR/RSC_ABCsims_20250722.R ${SLURM_ARRAY_TASK_ID} 100 ${USER}
