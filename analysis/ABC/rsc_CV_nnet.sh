#Neural network cross validations (10 reps each, 10 models = 100 total jobs) - 0.025 tolerance
cd /mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SHELL

cat modList.txt | while read MOD; do
echo '#!/bin/bash --login
##SBATCH -C [intel16|intel18|amd20]
#SBATCH -N 1 -c 1
#SBATCH -t 54:00:00
#SBATCH --mem 100G
#SBATCH -J rsc_CV_nnet

newgrp - Scribner_Lab

#module load  GCC/11.2.0  OpenMPI/4.1.1 R/4.1.2
module purge
module load GCC/12.2.0 OpenMPI/4.1.4-GCC-12.2.0 R/4.2.2-foss-2022b ICU/72.1-GCCcore-12.2.0

R_DIR="/mnt/ufs18/rs-012/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/"

source $R_DIR/setup_R.sh


cd /mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SHELL

Rscript rsc_CV_subset.R neuralnet '${MOD}' 5 0.025 ${SLURM_ARRAY_TASK_ID}

scontrol show job ${SLURM_JOB_ID}' > rsc_100k_CV_nnet_$MOD.sh 

sbatch -a 1-20 -o "/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/QSTAT/rsc_100k_CV_nnet_$MOD-%a.out" rsc_100k_CV_nnet_$MOD.sh 

done
