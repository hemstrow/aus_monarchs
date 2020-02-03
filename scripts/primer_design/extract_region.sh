#!/bin/bash -l
#SBATCH --mem=24G
#SBATCH -t 6-12:00:00
#SBATCH -J primer3

genome=~/genomes/Dapl_Zhan_v3_HiC.RN.fasta
chr=chr2
st=7820867
en=7881485
output=target_region.fasta

samtools faidx $genome ${chr}:${st}-${en} > $output
