#!/bin/bash/
#$ -j y
#$ -S /bin/bash 
#$ -V
#$ -l hostname=beast
#$ -pe multi 16
#$ -m e
#$ -M ly36@st-andrews.ac.uk

conda activate phastENV

#Use neutral/non-conserved model to estimate conserved 
cat chrom_names.txt | while read line;do
phastCons --msa-format MAF --target-coverage 0.25 --expected-length 12 --rho 0.4 --most-conserved $line.most-conserved_elements.gff CHROM_MAFS_PHAST/$line.forphast.maf nonconserved-4d.mod
done
