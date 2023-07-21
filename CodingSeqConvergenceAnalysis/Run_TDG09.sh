#!/bin/bash/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -l hostname=beast
#$ -pe multi 16


#EXAMPLE of script used to run analyses with TDG09. 



for phylip in /storage/home/users/ly36/FISH_WGA/Genomes/WGA_CYPRINODON/MafFilter/AA_alignments_TDG/*.phylip;do
java -cp tdg09/dist/tdg09.jar tdg09.Analyse -alignment $phylip -tree cyprinodon_tree.nw -groups OV VI -threads 16
done
