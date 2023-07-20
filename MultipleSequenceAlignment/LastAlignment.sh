#!/bin/bash/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -l hostname=slayer
#$ -pe multi 6
#$ -m e
#$ -M ly36@st-andrews.ac.uk

#An example of alignment of genomes to the reference genome performed for all 21 species in the analysis. 

#activate whole-genome alignment file.
conda activate WGA_env

#Align cyprinid fish genomes to g. multiradiatus genome using last aligner with three commands: lastdb, last-train, lastal.


# 1. Make a database for the reference genome, in this case we used G. multiradiatus. 
lastdb -P0 -uNEAR -R01 GM-assembly-lastdb GMassembly.fasta

#2. last-train is then used to improve alignment in the next step. last-train infers subsituitions, deletion and gap-rate and uses inferred parameters to aid alingment.
last-train -P0 --revsym --matsym --gapsym -E0.05 -C2 GM-assembly-lastdb Fheteroclitus.fas > GM_Fhetero.mat

#3. lastal is the alignment step where the target genome is aligned to the reference. Output is in maf format.
lastal -m50 -E0.05 -C2 -p GM_Fhetero.mat GM-assembly-lastdb Fheteroclitus.fas | last-split -m1 > GM_Fhetero.maf


