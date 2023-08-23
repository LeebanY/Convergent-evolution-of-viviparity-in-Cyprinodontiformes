## Non-coding analyses for Cyprinid convergent evolution of viviparity.

Here are scripts that were used to extract non-coding sequences, infer conserved non-coding sequences by first estimating the subsitution rate in four-fold degenerate sites across the 21 species, and then finally scripts to identify non-coding sequences that experiences accelerated sequence evolution in branches where the transition to viviparity occured. You will need install "WGAbed", "bedtools" and "Phast" to run the following scripts. 

1. MafToBed.sh and Gff_to_chrgff.sh are scripts to convert maf sequences into bed format with coordinates for all sequences, and to convert gff annotation file to chr specific gff files.
   
2. extract_4d_sites.sh : In this script, I am using gff and bed file specific chrs generated in previous scripts to extract four-fold degenerate sites using the msa-view from the PHAST program. The script then convert chr specific 4-fold-degenerate site files into ss format.
   
3. estimatenonconservedregions.sh : This script is used to output the most conserved elements in the genome -- this includes coding and non-coding regions. To do this, we use phastCons. Phastcons estimates most conserved regions of the genome using the REV model of substituition rates of 4-fold degenerate sites by asking which regions of the genome are involving significantly slower than 4-fold sites.

4. ExtractNonExonicElement.sh: Use output of estimatenonconservedregions.sh to extract only conserved non-coding elements from gff file and maf file.

5. Phyloacc.sh : A script to detect accelerated sequence divergence at branches where the transition to viviparity occured. Specifically tests for acceleration of sequence evolution at foreground branches rather than increased conservation, which is also possible with phast.

6. Noncoding_Regions_rapid_evolution_CNEE.R : Downstream analysis in R to perform exploratory analysis and to create figures.
   
 
