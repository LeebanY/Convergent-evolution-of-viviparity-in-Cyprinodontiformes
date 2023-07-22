#!/bin/bash/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -l hostname=beast
#$ -pe multi 8
#$ -m e
#$ -M ly36@st-andrews.ac.uk


#EXAMPLE -- A script to detect accelerated sequence divergence at branches where the transition to viviparity occured. 

# PhyloP used with a likelihood ratio test. --mode specifies that we want to test accelerated sequence divergence rather than conservation.


#conda deactivate
conda activate phastENV

cat chrom_names.txt | while read line;do
phyloP --msa-format MAF --method LRT --branch Gambusia_holbrooki-Poecilia_formosa,Girardinichthys_multiradiatus-Xenotaenia_resolanae --mode ACC --features GFFS_CNEE/$line.gff non_conserved_ancestors_annotated-4d.mod CHROM_MAFS_PHAST/$line.forphast.maf > CNEE_ACC/$line.phyloacc_NONEXONIC_results.txt
done
