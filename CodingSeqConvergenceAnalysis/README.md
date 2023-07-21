# Convergence analysis on coding sequence.

## TDG09

TDG09 is software developed to detect changes in selective constraint for different groups of lineage in a tree. For each amino acid in a multiple sequence alignment, it compares a null model where amino acid frequencies across a tree are homogenous, to a model with heterogenous amino acid frequencies between defined groups of lineages.

**CITATION: Tamuri AU, dos Reis M, Hay AJ, Goldstein RA (2009) Identifying Changes in Selective Constraints: Host Shifts in Influenza. PLoS Comput Biol 5(11): e1000564. doi:10.1371/journal.pcbi.1000564**
1. Run_TDG09.sh is an example script which illustrates how to run TDG09 (once installed), alongside an example phylip file and a required newick tree.
2. ProcessingTDGOutputs.sh is an example of commands used to process the outputs of the TDG09 analysis once finished.

## RERconverge.
RERconverge is a method that looks at rate of protein evolution across genes by comparing branch lengths of proteins.
**CITATION: https://github.com/nclark-lab/RERconverge and https://academic.oup.com/bioinformatics/article/35/22/4815/5514536**
1. /RERConvergeAnalysis.R : details how exactly the RER analysis was done -- follow the vignette at https://github.com/nclark-lab/RERconverge for specific details on how to run RERconverge for your dataset.
