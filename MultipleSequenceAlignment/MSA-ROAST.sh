#!/bin/bash/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -l hostname=phylo
#$ -pe multi 16
#$ -m e
#$ -M ly36@st-andrews.ac.uk

conda activate WGA_env

#EXAMPLE of script to perform multiple sequence alignment once netting and chaining has been performed.
 
maf_source=FISH_WGA/Genomes/CHAINING/FINAL_MAFS_22_2020/
OUT_destination=~/FISH_WGA/Genomes/WGA_CYPRINODON/ROASTED_CYPRINODON_ATTEMPT_3.maf
ref=Girardinichthys_multiradiatus

roast E=$ref "((((((Gambusia_holbrooki Gambusia_affinis) ((Xiphophorus_couchianus Xiphophorus_maculatus) Xiphophorus_hellerii)) ((Poeciliopsis_occidentalis Poeciliopsis_turrubarensis) Poeciliopsis_retropinna)) (((Poecilia_formosa Poecilia_mexicana) Poecilia_latipinna) Poecilia_reticulata)) Orestias_ascotanensis) (((Cyprinodon_nevadensis Cyprinodon_variegatus) Fundulus_heteroclitus) (((Girardinichthys_multiradiatus (Xenoophorus_captivus Goodea_atripinnis)) Xenotaenia_resolanae) Crenichthys_baileyi)))" $maf_source $OUT_destination
