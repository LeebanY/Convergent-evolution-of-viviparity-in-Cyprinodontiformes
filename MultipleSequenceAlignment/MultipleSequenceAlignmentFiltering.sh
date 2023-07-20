#!/bin/bash/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -l hostname=beast
#$ -pe multi 16

conda activate prankENV

#EXAMPLE script: Aligning and filtering msa for genes.

#Use masce v2 to implement a set of filtering guidelines.
#use mafft to do initial alignment because of speed.
mafft $input > $output
java -jar macse_v2.01.jar -prog trimNonHomologousFragments -seq $fas -out_NT $fas.NT.fas
#align sequences but balance alignment with how speed -- this can be tinkered with to ensure balance between good alignment and speed.
java -jar macse_v2.01.jar -prog alignSequences -seq $fas.NT.fas -max_refine_iter 0 -local_realign_init 0.05 -local_realign_dec 0.05 -out_NT $fas.alignedmacse.fas
#replace ! which are not supported into 'N'.
sed -i 's/!/N/g' $fas
#translate filtered amino acid sequences into nucleotide sequence. 
java -jar macse_v2.01.jar -prog translateNT2AA -seq $fas -out_AA $fas.AA.fasta
#use pal2nal to convert to multiple sequence alignment which is definitely inframe and remove mismatches and codons with gaps.
pal2nal.pl $fas.AA.fasta $fas -nomismatch -nogap -output fasta > $fas.pal2nal.fasta
