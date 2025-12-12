#!/bin/bash --login
#SBATCH --time=23:59:59
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=16G
#SBATCH --job-name=rsc_RF_pred
##SBATCH --error=../QSTAT/rsc_RF_pred.%J.stderr
#SBATCH --output=../QSTAT/rsc_RF_pred.%J.stdout

# Usage:  sbatch ../SHELL/rsc_RF_pred_subset.sh
# run from /mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/

module purge
module load GCC/12.2.0 OpenMPI/4.1.4-GCC-12.2.0 R/4.2.2-foss-2022b ICU/72.1-GCCcore-12.2.0

R_DIR="/mnt/ufs18/rs-012/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/"
OUT_DIR="/mnt/ufs18/rs-012/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/OUT2/"

source $R_DIR/setup_R.sh

cd $OUT_DIR

# Run model selection
$R_DIR/SRC/R-4.1.2-install/bin/Rscript ../SHELL/rsc_RF_pred_subset.R