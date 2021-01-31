#!/bin/bash/
#$ -j y
#$ -S /bin/bash 
#$ -V
#$ -l hostname=phylo
#$ -pe multi 1
#$ -m e
#$ -M ly36@st-andrews.ac.uk

#Extract coding regions from gff file. This is required to get 4d sites.

cat chrom_names.txt | while read line;do
less GM.augustus.filter_CDS.gff3 | grep "$line" > GFFS_coding/$line.gff
done
