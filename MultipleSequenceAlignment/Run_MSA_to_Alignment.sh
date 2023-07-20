#!/bin/bash/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -l hostname=slayer
#$ -pe multi 6
#$ -m e
#$ -M ly36@st-andrews.ac.uk

#Use maf filter to extract genes from multiple sequence alignment.

./maffilter param=MSAalignment_to_PerGeneAlignment.bpp


# Maf filter citation: Dutheil, J.Y., Gaillard, S. and Stukenbrock, E.H., 2014. MafFilter: a highly flexible and extensible multiple genome alignment files processor. BMC genomics, 15, pp.1-10.
