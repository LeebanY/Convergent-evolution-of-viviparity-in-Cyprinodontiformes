#!/bin/bash/
#$ -j y
#$ -S /bin/bash 
#$ -V
#$ -l hostname=phylo
#$ -pe multi 1
#$ -m e
#$ -M ly36@st-andrews.ac.uk

conda activate python2
module load bcftools


#Script used to split whole-genome alignment MAF file up into smaller MAF files by scaffold name. 

gzip ROASTED_CYPRINODON_ATTEMPT_singlecov_HEADTRIMMED.maf

cat chrom_names.txt | while read line;do
zgrep -v ^6 ROASTED_CYPRINODON_ATTEMPT_singlecov_HEADTRIMMED.maf.gz | WGAbed/maf_extract_ref_chr.py -c $line > CHROM_MAFS_PHAST/$line.forphast.maf
done
