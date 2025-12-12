#!/bin/bash --login
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --job-name=PE
##SBATCH --error=../QSTAT/pe.%J.stderr
#SBATCH --output=../QSTAT/pe.%J.stdout

MOD=$1
METH=$2
TOL=$3

# Usage: Rscript rsc_PE_subset.R SHD_bridge_pROW neuralnet 0.05

module purge
#module load  GCC/11.2.0 OpenMPI/4.1.1 R/4.1.2
module load GCC/12.2.0 OpenMPI/4.1.4-GCC-12.2.0 R/4.2.2-foss-2022b ICU/72.1-GCCcore-12.2.0

R_DIR="/mnt/ufs18/rs-012/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/"

source $R_DIR/setup_R.sh


$R_DIR/SRC/R-4.1.2-install/bin/Rscript rsc_PE_subset.R $MOD $METH $TOL
