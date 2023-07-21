#!/bin/bash/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -l hostname=beast
#$ -pe multi 16


#EXAMPLE of script used to run analyses with TDG09. 
#PREQUISITIES:
# TDG09 compares a null-model where amino acid frequencies at a given site are expected to be homogenous across branches in a tree, to a model where they are heterogenous with respect to defined groups of lineages. In this case, we specify two groups: viviparous lineages and oviparous lineages.
# Specifying branches requires producing MSA in phylip format where species contain a label (for example OV or VI, see example phylip file in this folder). 
# The labels included in the phylip file should also be replicated in the tree provided to TDG09 (again, see ExampleTDG09Tree.nwk for an example). 
# The tree, phylip file as well as group labels/names should then be provided to TDG09. 

for phylip in /storage/home/users/ly36/FISH_WGA/Genomes/WGA_CYPRINODON/MafFilter/AA_alignments_TDG/*.phylip;do
java -cp tdg09/dist/tdg09.jar tdg09.Analyse -alignment $phylip -tree cyprinodon_tree.nw -groups OV VI -threads 16
done


# CITATION: Tamuri AU, dos Reis M, Hay AJ, Goldstein RA (2009) Identifying Changes in Selective Constraints: Host Shifts in Influenza. PLoS Comput Biol 5(11): e1000564. doi:10.1371/journal.pcbi.1000564
