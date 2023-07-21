# Multiple sequence alignment for cyprinids. 
A set of scripts used in Yusuf et al (2022) to get from Fasta genome files to filtered NT and AA sequences for convergence analysis.

1. Last alignment (LastAlignment.sh) script to perform pairwise alignment of query sequences (cyprinid genomes) to the reference genome. See script for details on usage.
2. Chaining and netting (Chaining.sh): script to perform chaining and netting. Chains are sequences of blocks with matching query and reference sequence that are scored based on gappiness. Netting is the process of producing a hierarchical collection of chains which are ranked according to their gappiness. The best chains are contain the longest blocks with least gaps which may subsequently be filled in by chains with lower scores. See details in script. For more information see resource: http://genomewiki.ucsc.edu/index.php/Chains_Nets.
3. Perform sequence multiple genome alignment (ROAST-MSA.sh)
4. MSAalignment_to_PerGeneAlignment.bpp and Run_MSA_to_Alignment.sh: Extract protein-coding genes from the multiple sequence alignment.
5. MultipleSequenceAlignmentFiltering.sh: A script with utilising a combination of MAFFT (for fast alignment), MACSE (for codon-aware alignment with parameters that optimise speed) and PAL2NAL to convert sequences to nucleotide sequence.
6. ProteinFiltering.sh: A script used to filter MACSE AA outputs for convergence analysis.
