#!/bin/bash/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -l hostname=slayer
#$ -pe multi 16

#EXAMPLE script to perform branch tests for positive selection. Requires an coding sequence alignment, a tree with annotated branches and specifying which of those branches to test (Foreground).
# In this case the foreground branches are branches where the transition to viviparity occured.

conda activate prankENV

location=/storage/home/users/ly36/FISH_WGA/Genomes/WGA_CYPRINODON/MafFilter/Alignments_2/Convergent_genes_PAL2NAL

for fas in $location/*.fasta;do
hyphy absrel --alignment $fas -tree /storage/home/users/ly36/FISH_WGA/Genomes/WGA_CYPRINODON/MafFilter/HYPHY/annotated_cyprinodon_tree.txt --branches Foreground --CPU 16
done
