#!/bin/bash -l
#SBATCH --array=1-100
#SBATCH --mem=24G
#SBATCH -t 60-12:00:00
#SBATCH -J monRFrefine

# runs rf_refinement.R to run a refinement of a random forest model multiple times.

cd ~/monarch/aus_monarchs/aus_monarchs/results/rf_refinement
mkdir r_${SLURM_ARRAY_TASK_ID}
cd r_${SLURM_ARRAY_TASK_ID}

Rscript ~/monarch/aus_monarchs/aus_monarchs/scripts/run_rf_refinement.R ${SLURM_ARRAY_TASK_ID} rf_refinement

