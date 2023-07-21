#!/bin/bash/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -l hostname=beast
#$ -pe multi 16
#$ -m e
#$ -M ly36@st-andrews.ac.uk


#EXAMPLE SCRIPT -- In this script, I am using gff and bed file specific chrs generated in previous scripts to extract four-fold degenerate sites using the msa-view from the PHAST program.
# I then convert chr specific 4-fold-degenerate site files into ss format.
# these chr specific files are then aggregated to form one large file containing all extracted four fold degenerate sites and specifying the species in the MSA from which that data was originally extracted.
# Finally, a REV model is fit to model substituition rates of four-fold degenerate sites across the Cyprinodontiformes.

conda activate phastENV

cat chrom_names.txt | while read line;do
msa_view CHROM_MAFS_PHAST/$line.forphast.maf --in-format MAF --4d --out-format SS --4d --features GFFS_coding/$line.gff > 4d-codons.$line.ss
msa_view 4d-codons.$line.ss --in-format SS --out-format SS --tuple-size 1 > 4d-sites.$line.ss
msa_view --aggregate Girardinichthys_multiradiatus,Goodea_atripinnis,Fundulus_heteroclitus,Cyprinodon_nevadensis,Cyprinodon_variegatus,Xenotaenia_resolanae,Xenoophorus_captivus,Crenichthys_baileyi,Gambusia_holbrooki,Gambusia_affinis,Xiphophorus_couchianus,Xiphophorus_maculatus,Xiphophorus_hellerii,Poeciliopsis_turrubarensis,Poeciliopsis_retropinna,Poecilia_formosa,Poecilia_mexicana,Poecilia_latipinna,Poecilia_reticulata,Orestias_ascotanensis,Poeciliopsis_occidentalis 4d-sites.$line.ss > all-4d.sites.ss
phyloFit --tree "((((((Gambusia_holbrooki,Gambusia_affinis),((Xiphophorus_couchianus,Xiphophorus_maculatus),Xiphophorus_hellerii)),((Poeciliopsis_occidentalis,Poeciliopsis_turrubarensis),Poeciliopsis_retropinna)),(((Poecilia_formosa,Poecilia_mexicana),Poecilia_latipinna),Poecilia_reticulata)),Orestias_ascotanensis),(((Cyprinodon_nevadensis,Cyprinodon_variegatus),Fundulus_heteroclitus),(((Girardinichthys_multiradiatus,(Xenoophorus_captivus,Goodea_atripinnis)),Xenotaenia_resolanae),Crenichthys_baileyi)))" --msa-format SS --out-root nonconserved-4d all-4d.sites.ss
done
