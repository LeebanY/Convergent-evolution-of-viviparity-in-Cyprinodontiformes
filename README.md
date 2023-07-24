# Genomic signatures associated with transitions to viviparity in Cyprinodontiformes.
Comparative genomics project examining the role of convergent evolution in the independent origin of viviparity in Cyprinodontiformes. 

## Aims of the analysis.
Viviparity has evolved multiple times, but the evolution of viviparity, and importantly, whether the same genes are always recruited is still unknown. Here, we assess in two clades of fish that have evolved viviparity, whether the same genes have convergently evolved. 

## Brief description of the scripts.
This GitHub repo contains scripts necessary to reperform the analysis. Steps include, aligning, chaining, netting, multiple sequence alignment of public and de-novo genomes, extraction of different genomic features, alignment of orthologs, analysis of convergent evolution at amino acid sites, inferring phylogenetic trees for every sequence, performing analysis to examine evolutionary rates for every protein and finally examining convergent evolution in non-coding regions. 

## Order of analysis.
1. Multiple Sequence Alignment (MultipleSequenceAlignment).
     - This contains scripts required to perform multiple sequence alignment of genomes. Requires KentUtils, TBA and MULTIZ programs to function.
2. Protein-coding sequence test for convergent evolution (CodingSeqConvergent).
     - This includes scripts that filter and align protein-coding genes and scripts to perform convergence analysis. Requires MAFFT, MASCE, MafFilter, Divvier, TDG09, Phangorn and RERconverge.
3. Non-coding sequence test for convergent evolution (NoncodingSequenceAnalysis).
      - This contains scripts to extracting relevant non-coding loci and perform non-coding analysis. Requires PHAST and WGAbed.
4. Checks for molecular convergence and ILS (MoleConvergenceILSTests)






