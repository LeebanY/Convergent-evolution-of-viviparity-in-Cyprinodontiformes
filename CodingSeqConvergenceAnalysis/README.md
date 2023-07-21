# Convergence analysis on coding sequence.

## TDG09

TDG09 is software developed to detect changes in selective constraint for different groups of lineage in a tree. For each amino acid in a multiple sequence alignment, it compares a null model (H0) where amino acid frequencies across a tree are homogenous, to a model with heterogenous amino acid frequencies between defined groups of lineages (H1).

**CITATION: Tamuri AU, dos Reis M, Hay AJ, Goldstein RA (2009) Identifying Changes in Selective Constraints: Host Shifts in Influenza. PLoS Comput Biol 5(11): e1000564. doi:10.1371/journal.pcbi.1000564**
1. Run_TDG09.sh is an example script which illustrates how to run TDG09 (once installed), alongside an example phylip file and a required newick tree.
       -  Prerequisites include making sure to include group-specific labels to the phylip and newick tree used for analysis. 
3. ProcessingTDGOutputs.sh is an example of commands used to process the outputs of the TDG09 analysis once finished. The end file will be a list Amino acid alignments followed by any sites with significant p-values following a likelihood ratio test between H0 and H1. 

## RERconverge.
RERconverge is a method that looks at rate of protein evolution across genes by comparing branch lengths of proteins.
**CITATION: https://github.com/nclark-lab/RERconverge and https://academic.oup.com/bioinformatics/article/35/22/4815/5514536**
1. /RERConvergeAnalysis.R : details how exactly the RER analysis was done -- follow the vignette at https://github.com/nclark-lab/RERconverge for specific details on how to run RERconverge for your dataset.
      - Rerconverge first requires trees to produced with the same tree topology but where branch lengths are estimated based on each amino acid multiple-seq alignment. This is done using phangorn with LG as an AA subsituition model.
      - Trait of interest (viviparity) is labelled in the tree by specifying (1) at the foreground branches of interest.

RERconverge outputs a correlation between the relative evolutionary rate and the trait of interest.
      - 
