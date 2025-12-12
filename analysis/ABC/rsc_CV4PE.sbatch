#!/bin/bash --login
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --job-name=CV4PE
##SBATCH --error=../QSTAT/cv4pe.%J.stderr
#SBATCH --output=../QSTAT/cv4pe.%J.stdout

MOD=$1
METH=$2
TOL=$3


module purge
#module load  GCC/11.2.0 OpenMPI/4.1.1 R/4.1.2
module load GCC/12.2.0 OpenMPI/4.1.4-GCC-12.2.0 R/4.2.2-foss-2022b ICU/72.1-GCCcore-12.2.0

R_DIR="/mnt/ufs18/rs-012/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/"

source $R_DIR/setup_R.sh


#Rscript rsc_CV4PE.R Grp2_bridge_indGrp1-5 neuralnet 10 0.05

Rscript rsc_CV4PE.R $MOD $METH 10 $TOL
