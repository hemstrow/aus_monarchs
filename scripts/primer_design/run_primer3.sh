#!/bin/bash -l
#SBATCH --mem=24G
#SBATCH -t 6-12:00:00
#SBATCH -J primer3

parms=parms
fasta=target_region.fasta
outfile=MonarchMigRegion

Rscript get_primer_intervals.R

Rscript prepare_primer3_parms.R $parms primer_regions.txt $fasta $outfile

primer3_core --format_output --strict_tags --output=${oufile}_primer3_summary.txt < parms_${oufile}

