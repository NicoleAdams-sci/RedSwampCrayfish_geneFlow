#!/bin/bash --login
#SBATCH --time=00:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --job-name=mnlog
##SBATCH --error=mnlog.%J.stderr
#SBATCH --output=../QSTAT/mnlog.%J.stdout

TOL=$1

module purge
module load GCC/12.2.0 OpenMPI/4.1.4-GCC-12.2.0 R/4.2.2-foss-2022b ICU/72.1-GCCcore-12.2.0

R_DIR="/mnt/ufs18/rs-012/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/"

source $R_DIR/setup_R.sh

$R_DIR/SRC/R-4.1.2-install/bin/Rscript $R_DIR/SHELL/rsc_modsel_subset_mnlog.R $TOL
