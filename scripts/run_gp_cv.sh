#!/bin/bash -l
#SBATCH --array=1-203
#SBATCH --mem=20G
#SBATCH -t 6-12:00:00
#SBATCH -J monGPcv

# runs gp_cv.R to do a leave one out cv for every (clean) sample in genotypes.geno

cd ~/monarch/aus_monarchs/aus_monarchs/results/gp_leave_one_out
mkdir r_${SLURM_ARRAY_TASK_ID}
cd r_${SLURM_ARRAY_TASK_ID}

Rscript ~/monarch/aus_monarchs/aus_monarchs/scripts/gp_cv.R ${SLURM_ARRAY_TASK_ID} ../BayesB_leave_one_out_cv.txt

