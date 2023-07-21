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

#EXAMPLE script -- convert maf MSA to a chr-specific bed file for each maf to use later on to extract scaffold/chromosome specific sequences.

gzip ROASTED_CYPRINODON_ATTEMPT_singlecov_HEADTRIMMED.maf

cat chrom_names.txt | while read line;do
#python WGAbed/maf_to_bed.py -i ROASTED_CYPRINODON_ATTEMPT_singlecov_HEADTRIMMED.maf.gz -r Girardinichthys_multiradiatus -c $line | sort -k1,1 -k2,2n | bgzip -c > Chrom_Mafs_data/$line.wga.bed.gz
zgrep -v ^6 ROASTED_CYPRINODON_ATTEMPT_singlecov_HEADTRIMMED.maf.gz | WGAbed/maf_extract_ref_chr.py -c $line > CHROM_MAFS_PHAST/$line.forphast.maf
done
