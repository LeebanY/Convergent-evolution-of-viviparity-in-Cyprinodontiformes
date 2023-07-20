#!/bin/bash/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -l hostname=phylo
#$ -pe multi 16

conda activate prankENV

#Filter protein alignments using divvier partial. Ensuring that columns or amino acids that are kept are columns with no gaps (i.e. no missing data in any species).

divvier partial -divvygap -mincol 21 $fas
