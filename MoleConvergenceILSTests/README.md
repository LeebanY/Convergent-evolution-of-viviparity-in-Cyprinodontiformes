## Analysis of ILS for molecular convergence in RER. 
RER relies on a fixed topology, where branch length differences across different sequences/genes at particular branches are assessed. To test whether RER is sensitive to ILS, we inferred alternative species tree topologies using IQTREE that differed particularly at foreground branches (FocalBranchGoodeid1TRIAT.nwk, for example). We then repeated our RERconverge analysis using these alternative topologies and assessed the overlap of genes showing convergent sequence divergence across different topologies.

RER_shell.sh was used to run RER_revision_analysis.R on a SLURM cluster. Contains details on how to infer trees using phangorn using these alternate topologies.

