#!/bin/bash -l
#SBATCH -t 1-12:00:00
#SBATCH --mem 16G
#SBATCH -c 16

bamlist=../aus_monarchs/data/gbamlist_minus_bads.txt
out=association_out_likely
ybin=../aus_monarchs/data/likely_ybin.txt

minInd=102

~/angsd/angsd -P 16 -bam $bamlist -out $out -minQ 20 -minMapQ 20 -GL 1 -SNP_pval 1e-12 -doMajorMinor 1 -doMaf 2 -minMaf 0.05 -minInd $minInd -doAsso 1 -yBin ${ybin}

