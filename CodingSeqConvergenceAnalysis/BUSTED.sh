#!/bin/bash/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -l hostname=slayer
#$ -pe multi 16

#EXAMPLE script for performing branch-site tests asking if there is evidence for positive selection at any site across any branch in the tree.

conda activate prankENV

location=/storage/home/users/ly36/FISH_WGA/Genomes/WGA_CYPRINODON/MafFilter/Alignments_2/Convergent_genes_PAL2NAL

for fas in $location/*.fasta;do
hyphy busted --alignment $fas -tree /storage/home/users/ly36/FISH_WGA/Genomes/WGA_CYPRINODON/MafFilter/HYPHY/annotated_cyprinodon_tree.txt --CPU 16
done
