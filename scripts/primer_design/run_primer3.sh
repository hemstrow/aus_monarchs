#!/bin/bash -l
#SBATCH --mem=24G
#SBATCH -t 6-12:00:00
#SBATCH -J primer3

primer3_core --format_output --strict_tags --output=primer3_output.txt < parms

