#!/bin/bash/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -l hostname=phylo
#$ -pe multi 1
#$ -m e
#$ -M ly36@st-andrews.ac.uk


cat chrom_names.txt | while read line;do
less All_Scaffold.NONEXONICELEMENTS.most-conserved_elements.phast.gff | grep "$line" > GFFS_CNEE/$line.gff
done
