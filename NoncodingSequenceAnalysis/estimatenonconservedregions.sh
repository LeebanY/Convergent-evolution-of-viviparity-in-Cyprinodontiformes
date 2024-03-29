#!/bin/bash/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -l hostname=beast
#$ -pe multi 16
#$ -m e
#$ -M ly36@st-andrews.ac.uk

#EXAMPLE script -- This script is used to output the most conserved elements in the genome -- this includes coding and non-coding regions. To do this, we use phastCons. 
#phastcons estimates most conserved regions of the genome using the REV model of substituition rates of 4-fold degenerate sites by asking which regions of the genome are involving significantly slower than 4-fold sites. These are then outputted.


conda activate phastENV

cat chrom_names.txt | while read line;do
phastCons --msa-format MAF --target-coverage 0.25 --expected-length 12 --rho 0.4 --most-conserved $line.most-conserved_elements.gff CHROM_MAFS_PHAST/$line.forphast.maf nonconserved-4d.mod
done
