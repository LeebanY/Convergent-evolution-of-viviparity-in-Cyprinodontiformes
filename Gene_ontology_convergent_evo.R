library(tidyverse)
library(viridis)
library(phangorn)
library(data.table)
library(gt)
library(ggrepel)
library(ggpubr)
library(stringr)
#The following dataset was produced by taking the diamond blast results, extracting the accession number
#Then, accession numbers were converted into gene names via uniprot and duplicates were removed.
#Then, gene names were put into PANTHER (web) for overepresentation analysis with zebrafish as REF
convergent_gene_set_Drerio <- read.table("Gene_ontology_PANTHER_overrepresentation_danioref.txt",
                                  fill = TRUE , header = T)
convergent_gene_set_Hsap <- read.table("Gene_ontology_PANTHER_human_background_convergent.txt",
                                         fill = TRUE , header = T)


tiff("GO_convergent_geneset_against_Drerio.tiff", units="in", width=7, height=5, res=300)
ggplot(top_10_convergent_zebra, aes(x=Fold_enrichment, y=reorder(Gene_ontology_biological_process, Fold_enrichment), colour=FDR))+
  geom_point(size=top_10_convergent_zebra$Fold_enrichment*1.2)+scale_colour_viridis(option='D', discrete = F)+
  theme_minimal()+xlab('Fold Enrichment')+ylab('Gene ontology (Biological process - zebrafish background)')+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))
dev.off()

top_10_convergent_human<-top_n(convergent_gene_set_Hsap, 10, Fold_enrichment)
top_10_convergent_zebra<-top_n(convergent_gene_set_Drerio, 10, Fold_enrichment)
tiff("GO_convergent_geneset_against_Hsapiens.tiff", units="in", width=7, height=5, res=300)
ggplot(top_10_convergent_human, aes(x=Fold_enrichment, y=reorder(Gene_ontology_biological_process, Fold_enrichment), colour=FDR))+
  geom_point(size=top_10_convergent_human$Fold_enrichment*0.5)+scale_colour_viridis(option='D', discrete = F)+
  theme_minimal()+xlab('Fold Enrichment')+ylab('Gene ontology (Biological process - human background)')+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))
dev.off()




convergent_gene_set_Drerio$Gene_ontology_biological_process<-gsub("_", " ", convergent_gene_set_Drerio$Gene_ontology_biological_process)
top_30_convergent_human$Gene_ontology_biological_process<-gsub("_", " ", top_30_convergent_human$Gene_ontology_biological_process)
top_10_convergent_human$Gene_ontology_biological_process<-gsub("_", " ", top_10_convergent_human$Gene_ontology_biological_process)


tiff("gene_annotations_convergent_amino_acid_change_hsap_drer.tiff", units="in", width=11, height=5, res=300)
ggarrange(hsap, Drer, common.legend = T, font.label = list(size= 10, face='bold'), legend = "bottom" )
dev.off()



################ Trees for RERconverge #############
library(ape)
test_tree <- read.tree(file='cyprinodon_gene_tree_test.txt')
labelled_tree <- read.tree(file='Labelled_viviparity_cyp_tree.txt')
labelled_tree <- read.tree(file='Testing_seaview')
test_tree
unroot_test_tree<-root(test_tree, node=27)
plot(unroot_test_tree)+nodelabels(cex=0.3) # Both of them are the same.
plot(test_tree)+nodelabels(cex=0.3) #Same as unrooted
plot(labelled_tree)
plot()

##############**** RERCONVERGE ****#############
#RERconverge is a method that looks at rate of protein evolution across genes by comparing branch lengths of proteins.
#This is complimentary to looking for changes at a particular site. 
library(devtools)
install_github("nclark-lab/RERconverge", dependencies = T)
library(RERconverge)
#FINALLY RERCONVERGE HAS BEEN INSTALLED!!
rerpath = find.package('RERconverge')
#read in tree file
treefileconverge = "Phangorn_trees_AA" 
Treesconverge=readTrees(paste(treefileconverge), max.read = 17000)
########Produce your trees using phangorn
#estimatePhangornTreeAll(alndir = 'Fasta_AA_alignments', treefile = 'Cyprino_phangorn_tree.txt', output.file = 'Phangorn_trees_AA', format='fasta', type='AA')


read.aa('Fasta_AA_alignments/21Cyprinodontiformes10000fGirMul_m9_00083_085313538536967_AA.partial.fas.phylip.fasta',format='fasta')
##################
master<-Treesconverge$masterTree
spp_names<-master$tip.label
#Now get evolutionary rates that will be corrected for heteroscedasity.
fishviviparity= getAllResiduals(Treesconverge, useSpecies = spp_names,
                          transform = "sqrt", weighted = T, scale = T)#this works
#Save this just so you don't have to re-do.
saveRDS(fishviviparity, file="fishviviparityRER.rds") 
fishviviparity<-readRDS(fishviviparity, file='fishviviparityRER.rds')
# Just to test plot evolutionary rate for one gene.
#21Cyprinodontiformes18640fGirMul_m9_00005_066575736664312_AA.partial.fas.phylip  
outgroup <- c("OV_Oriestias_ascotanensis")
par(mfrow=c(1,2))
spp_names <-c("VI_Gambusia_holbrooki","VI_Xiphophorus_couchianus","VI_Xiphophorus_maculatus","VI_Xiphophorus_hellerii",
             "VI_Poeciliopsis_occidentalis","VI_Poeciliopsis_turrubarensis", "VI_Poeciliopsis_retropinna",
             "VI_Poecilia_formosa","VI_Poecilia_mexicana","VI_Poecilia_reticulata","VI_Poecilia_latipinna",
             "VI_Xenoophorus_captivus","VI_Girardinichthys_multiradiatus","VI_Xenotaenia_resolanae","VI_Goodea_atripinnis",
             "VI_Gambusia_affinis")
xtree='21Cyprinodontiformes18640fGirMul_m9_00005_066575736664312_AA.partial.fas.phylip '

gene_tree=plotTreeHighlightBranches(Treesconverge$trees$`21Cyprinodontiformes10000fGirMul_m9_00083_085313538536967_AA.partial.fas.phylip`, outgroup=outgroup,
                                  hlspecies=c("VI_Gambusia_holbrooki","VI_Xiphophorus_couchianus","VI_Xiphophorus_maculatus","VI_Xiphophorus_hellerii","VI_Poeciliopsis_occidentalis","VI_Poeciliopsis_turrubarensis", "VI_Poeciliopsis_retropinna", "VI_Poecilia_formosa","VI_Poecilia_mexicana","VI_Poecilia_reticulata","VI_Poecilia_latipinna","VI_Xenoophorus_captivus","VI_Girardinichthys_multiradiatus","VI_Xenotaenia_resolanae","VI_Goodea_atripinnis", "VI_Gambusia_affinis"), hlcols="red",
                                  main='Gene_Tree')
avgtree=plotTreeHighlightBranches(master, outgroup=outgroup,
                                  hlspecies=c("VI_Gambusia_holbrooki","VI_Xiphophorus_couchianus","VI_Xiphophorus_maculatus","VI_Xiphophorus_hellerii","VI_Poeciliopsis_occidentalis","VI_Poeciliopsis_turrubarensis", "VI_Poeciliopsis_retropinna", "VI_Poecilia_formosa","VI_Poecilia_mexicana","VI_Poecilia_reticulata","VI_Poecilia_latipinna","VI_Xenoophorus_captivus","VI_Girardinichthys_multiradiatus","VI_Xenotaenia_resolanae","VI_Goodea_atripinnis", "VI_Gambusia_affinis"), hlcols="red", 
                                  main="Average tree")
fishviviparity = readRDS("fishviviparityRER.rds")
plot(master)
#Binary trait analysis
trait_tree=read.tree("Cyprino_phangorn_tree_TRIATBINARY.txt")
phenotype_viviparity=tree2Paths(trait_tree, Treesconverge)
phenotype_viviparity_null1=tree2Paths(trait_tree, Treesconverge)
corViviparity=correlateWithBinaryPhenotype(fishviviparity, phenotype_viviparity, min.sp=21, min.pos=2)
Significant_correlation <- corViviparity %>% filter(P < 0.05) #532 significant correlations.
write.table(Significant_correlation, file = "RER_significant_correlations_genes.txt", sep='\t')
table(sign(Significant_correlation$Rho))
#Negative correlations -- what do they mean
gene_tree=plotTreeHighlightBranches(Treesconverge$trees$`21Cyprinodontiformes2503fGirMul_m9_00091_016201321672407_AA.partial.fas.phylip`, outgroup=outgroup,
                                    hlspecies=c("VI_Gambusia_holbrooki","VI_Xiphophorus_couchianus","VI_Xiphophorus_maculatus","VI_Xiphophorus_hellerii","VI_Poeciliopsis_occidentalis","VI_Poeciliopsis_turrubarensis", "VI_Poeciliopsis_retropinna", "VI_Poecilia_formosa","VI_Poecilia_mexicana","VI_Poecilia_reticulata","VI_Poecilia_latipinna","VI_Xenoophorus_captivus","VI_Girardinichthys_multiradiatus","VI_Xenotaenia_resolanae","VI_Goodea_atripinnis", "VI_Gambusia_affinis"), hlcols="red",
                                    main='Gene_Tree')
avgtree=plotTreeHighlightBranches(master, outgroup=outgroup,
                                  hlspecies=c("VI_Gambusia_holbrooki","VI_Xiphophorus_couchianus","VI_Xiphophorus_maculatus","VI_Xiphophorus_hellerii","VI_Poeciliopsis_occidentalis","VI_Poeciliopsis_turrubarensis", "VI_Poeciliopsis_retropinna", "VI_Poecilia_formosa","VI_Poecilia_mexicana","VI_Poecilia_reticulata","VI_Poecilia_latipinna","VI_Xenoophorus_captivus","VI_Girardinichthys_multiradiatus","VI_Xenotaenia_resolanae","VI_Goodea_atripinnis", "VI_Gambusia_affinis"), hlcols="red", 
                                  main="Average tree")


#assess distribution of P-values
hist(corViviparity$P, breaks=15, xlab="Kendall P-value", 
     main="P-values for correlation between 200 genes and marine environment")
#P-values look bimodal -- maybe this is concerning.

#Check if one species is driving the rate of evolutionary change
x=charpaths
y=mamRERw['TTN',]
pathnames=namePathsWSpecies(toyTrees$masterTree) 
names(y)=pathnames
plot(x,y, cex.axis=1, cex.lab=1, cex.main=1, xlab="Weight Change", 
     ylab="Evolutionary Rate", main="Gene ADSS Pearson Correlation",
     pch=19, cex=1, xlim=c(-2,2))
text(x,y, labels=names(y), pos=4)
abline(lm(y~x), col='red',lwd=3)
hist(corViviparity$P, breaks=15, xlab="Kendall P-value", 
     main="P-values for correlation between 200 genes and marine environment")

####### NULL/CONTROL model
Null_traittree_1=read.tree("Cyprino_phangorn_tree_TRIATBINARY_Null1.txt")
plot(Null_traittree_1)
NULL_viviparity=tree2Paths(Null_traittree_1, Treesconverge)
NULLCORViviparity=correlateWithBinaryPhenotype(fishviviparity, NULL_viviparity, min.sp=21, min.pos=2)
Significant_correlation_null1 <- NULLCORViviparity %>% filter(P < 0.05) 
table(sign(Significant_correlation_null1$Rho))
# 668 convergently evolving genes in empirical null model.

#Now run another null model - 3
Null_traittree_2=read.tree("Cyprino_phangorn_tree_TRIATBINARY_Null2.txt")
plot(Null_traittree_2)
NULL_viviparity_2=tree2Paths(Null_traittree_2, Treesconverge)
NULLCORViviparity_2=correlateWithBinaryPhenotype(fishviviparity, NULL_viviparity_2, min.sp=21, min.pos=2)
Significant_correlation_null2 <- NULLCORViviparity_2 %>% filter(P < 0.05) 
table(sign(Significant_correlation_null2$Rho))
#746 convergently evolving genes in empirical null model 2

#Now run another null model - 3
Null_traittree_3=read.tree("Cyprino_phangorn_tree_TRIATBINARY_Null3.txt")
plot(Null_traittree_3)
NULL_viviparity_3=tree2Paths(Null_traittree_3, Treesconverge)
NULLCORViviparity_3=correlateWithBinaryPhenotype(fishviviparity, NULL_viviparity_3, min.sp=21, min.pos=2)
Significant_correlation_null3 <- NULLCORViviparity_3 %>% filter(P < 0.05) 
table(sign(Significant_correlation_null3$Rho))

#574 convergently evolving genes in empirical null model 3.

##Mutate all datasets to include directionality.
Significant_correlation<-Significant_correlation %>%
  mutate(Direction=case_when(
    Rho < 0 ~ "Negative",
    Rho > 0 ~ "Positive"
  ))

Significant_correlation_null1<-Significant_correlation_null1 %>%
  mutate(Direction=case_when(
    Rho < 0 ~ "Negative",
    Rho > 0 ~ "Positive"
  ))

Significant_correlation_null2<-Significant_correlation_null2 %>%
  mutate(Direction=case_when(
    Rho < 0 ~ "Negative",
    Rho > 0 ~ "Positive"
  ))

Significant_correlation_null3<-Significant_correlation_null3 %>%
  mutate(Direction=case_when(
    Rho < 0 ~ "Negative",
    Rho > 0 ~ "Positive"
  ))



table(sign(Significant_correlation$Rho))
table(sign(Significant_correlation_null2$Rho))
table(sign(Significant_correlation_null3$Rho))

#Clearly no evidence for excess convergence at sequence level. 
tiff("Convergent_distribution_RER.tiff", units="in", width=9, height=6, res=300)
ggplot(Significant_correlation, aes(P, fill = Direction))+
  geom_histogram(binwidth = 0.002, alpha=0.8, col='black')+theme_pubclean()+
  labs(x='P-value (p<0.05)',
       y='Frequency')+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
dev.off()

Null_Experimental_RER_Corrs <- read.csv('Accelerated_Decelerated_Null_RER.csv')

tiff("Convergent_distribution_RER_NULLS_ACCDEC.tiff", units="in", width=4, height=6, res=300)
ggplot(Null_Experimental_RER_Corrs, aes(x=No_De_Genes, y=Correlations, fill=No_De_Genes))+
  geom_col(colour='black')+theme_pubr()+scale_fill_discrete(direction=-1)+
  facet_wrap(~Tests)+
  labs(x='Direction of sequence divergence', 
       y='Number of genes',
       fill='Direction of sequence divergence')+
  theme(axis.text=element_text(size=7.5), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position="none")
dev.off()


#Measure strength of divergence in null compared to foreground.

Significant_correlation<-Significant_correlation %>%
  mutate(Test='Experimental')
Significant_correlation_null1<-Significant_correlation_null1 %>%
  mutate(Test='Null 1')
Significant_correlation_null2<-Significant_correlation_null2 %>%
  mutate(Test='Null 2')
Significant_correlation_null3<-Significant_correlation_null3 %>%
  mutate(Test='Null 3')

RER_Experimental_Nulls_Altogether<-do.call("rbind", list(Significant_correlation, 
                      Significant_correlation_null1,
                      Significant_correlation_null2,
                      Significant_correlation_null3))

gene_con_lm<-lm(value ~ Chromosome, data = gene_con_melt)
gene_con_aov<-aov(value ~ Chromosome, data = gene_con_melt)

anova(gene_con_lm)
TukeyHSD(gene_con_aov)
RER_Experimental_null_model<-lm( Correlations ~ Tests,data = Null_Experimental_RER_Corrs)
anova(RER_Experimental_null_model)
RER_Exp_Null_aov <- aov(RER_Experimental_null_model)
TukeyHSD(RER_Exp_Null_aov)


tiff("Convergent_distribution_RER_NULLS_ACCDEC_violin.tiff", units="in", width=6, height=6, res=300)
ggplot(RER_Experimental_Nulls_Altogether, aes(x=Test, y=Rho, fill=Direction))+
  geom_violin()+stat_summary(fun = mean,
                             geom = "crossbar", 
                             width = 0.5)+
  labs(x='Tests conducted', 
       y='Relative evolutionary rate (Rho)',
       fill='')+geom_hline(yintercept = 0, linetype='dashed')+
  theme(legend.position="none")+theme_minimal()+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))

dev.off()
################ PLOTS, GET AVERAGE EVOLUTIONARY RATE FOR ACCELERATED AND DECELERATED GENES
fishviviparity_2 <- data.frame(fishviviparity)
fishviviparity_FILENAMES <- tibble::rownames_to_column(fishviviparity, "Filenames")


######## convert rownames to columns for significant
Significant_correlation_null$filenames <- rownames(Significant_correlation_null)
Significant_correlation$filenames <- rownames(Significant_correlation)
fishviviparity$filenames <- rownames(fishviviparity)




############### EXAMPLE FOR TALK -- GENES WITH ACCELERATED RER INVOLVED IN DEVELOPMENT ######
# tiff("SMAD.tiff", units="in", width=6, height=6, res=300)
# treePlotRers(treesObj=Treesconverge, rermat=fishviviparity, index="21Cyprinodontiformes16168fGirMul_m9_00029_01449655614503971_AA.partial.fas.phylip", 
#                             type="c",figwid=0.6) #SMAD2
# dev.off()
# tiff("NOG.tiff", units="in", width=6, height=6, res=300)
# treePlotRers(treesObj=Treesconverge, rermat=fishviviparity, index="21Cyprinodontiformes3148fGirMul_m9_00070_051101765111101_AA.partial.fas.phylip", 
#                             type="c",figwid=0.6) #NOG
# dev.off()
# 
# par(mar=c(1,1,1,1))
# #match
# Both_Significant_Correlations_Experimental_Null<-Significant_correlation[(Significant_correlation$filenames%in%Significant_correlation_null$filenames),]
# glimpse(Both_Significant_Correlations_Experimental_Null$filenames)
# 
# Positive_significant_correlations<-Significant_correlation %>%
#   filter(Rho > 0)
# 
# Negative_significant_correlations<-Significant_correlation %>%
#   filter(Rho < 0)
# 
# 
# 
# RER_genenames_morethan1<-RER_genenames %>%
#   filter(fpkm_mean > 1)
# 
# 

############################## BRANCH_SITE TESTS WITH RER CONVERGENCE ----#########

HYPHY_RER <- read.csv('Summary_of_ABSREL_BUSTED.csv')
HYPHY_RER <- na.omit(HYPHY_RER)

FDR_pvals_BUSTED<-p.adjust(HYPHY_RER$pvalue_BUSTED, method = 'BH', n=57)
HYPHY_RER<-cbind(HYPHY_RER,FDR_pvals_BUSTED)
FDR_pvals_RELAX<-p.adjust(HYPHY_RER$RELAX.P.value, method = 'BH', n=57)
HYPHY_RER<-cbind(HYPHY_RER,FDR_pvals_RELAX)


Busted_sig<-HYPHY_RER %>%
  filter(FDR_pvals_BUSTED < 0.05)

length(Busted_sig)


HYPHY_RER %>%
  filter(Found_in_either_Poeciliopsis=='Yes' & Found_in_Mammal_plac=='Yes' )

Interesting_genes<-HYPHY_RER %>%
  select(Genename, FDR_pvals_BUSTED, FDR_pvals_RELAX, Found_in_either_Poeciliopsis, Found_in_Mammal_plac, RELAX...Evidence.for.selection.) %>%
  filter(Found_in_either_Poeciliopsis=='Yes' & Found_in_Mammal_plac=='Yes' & FDR_pvals_BUSTED < 0.05 )

colnames(Interesting_genes) <- c("Genename", "P-value (BUSTED)", 'P-value (RELAX)', 'Found in maternal follicle of Poeciliopsis species','Found in mammalian placentas','Evidence of changes in selection pressure')
Interesting_genes[[1]] <- toupper(Interesting_genes[[1]])

Interesting_genes_gt<-gt(Interesting_genes)
gtsave(Interesting_genes_gt, '19_genes_Foundinplacenta_positiveselec_convergence.png')

neg_corrs_hyphy <- HYPHY_RER %>%
  filter(Correlation_coef_RER < 0)
pos_corrs_hyphy <- HYPHY_RER %>%
  filter(Correlation_coef_RER > 0)

ggplot(pos_corrs_hyphy, aes(BUSTED_proportion_div_sel, y=Correlation_coef_RER))+
  geom_point(aes(colour=pvalue_BUSTED))+theme_bw()+
  scale_colour_continuous()

ggplot(neg_corrs_hyphy, aes(BUSTED_proportion_div_sel, y=Correlation_coef_RER))+
  geom_point(aes(colour=pvalue_BUSTED))+theme_bw()+
  scale_colour_continuous()

ggplot(HYPHY_RER, aes(x=Genename, y=BUSTED_proportion_div_sel,
                      fill=Direction, size=pvalue_BUSTED))+
  geom_point()+theme_bw()+scale_fill_discrete()+
  coord_flip()


tiff("BUSTED_results_RERordered.tiff", units="in", width=6, height=8, res=300)
HYPHY_RER %>%
  arrange(Direction) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Genename=factor(Genename, levels=Genename)) %>%   # This trick update the factor levels
  ggplot( aes(x=Genename, y=BUSTED_proportion_div_sel, color=Direction, linetype=Significance)) +
  geom_segment( aes(xend=Genename, yend=0)) +
  geom_point(size=4) + ylab('Proportion of sites under diversifying selection (%)')+
  coord_flip() +
  theme_minimal() +
  xlab("")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
dev.off()

HYPHY

####### RER WITH GENENAMES PLOTTING #######
RER_genenames <- read.csv('RER_CORRELATIONS_WITH_GENENAMES.csv',header = T, fill=T)

Positive_significant_correlations<-RER_genenames %>%
  filter(Rho > 0)

Negative_significant_correlations<-RER_genenames %>%
  filter(Rho < 0)

ggplot(RER_genenames, aes(x=Rho, y=P))+
  geom_point()

Positive_significant_correlations$fpkm_mean <- Summary_rer_nodups$fpkm_mean[ match(Positive_significant_correlations$Genenames, Summary_rer_nodups$human_alignment_symbol)]
Negative_significant_correlations$fpkm_mean <- Summary_rer_nodups$fpkm_mean[ match(Negative_significant_correlations$Genenames, Summary_rer_nodups$human_alignment_symbol)]
RER_genenames$fpkm_mean <- Summary_rer_nodups$fpkm_mean[ match(RER_genenames$Genenames, Summary_rer_nodups$human_alignment_symbol)]
Positive_significant_correlations <- na.omit(Positive_significant_correlations)
Negative_significant_correlations <- na.omit(Negative_significant_correlations)
RER_genenames<-na.omit(RER_genenames)
summary(RER_genenames)

#Take all the RER genes with annotations - accelerated and conserved
Randomly_sampled_genes_comparison_acceleratedRER <- Background %>%
  slice_sample(n=296)

Replicated_Random_FPKM_rer<-replicate(1000, {Background %>%
  slice_sample(n=296) %>%
  summarise(mean(fpkm_mean))})

summary(Replicated_Random_FPKM_rer)

Replicated_Random_FPKM_rer<- data.frame(matrix(unlist(Replicated_Random_FPKM_rer), nrow=length(Replicated_Random_FPKM_rer), byrow=T))
names(Replicated_Random_FPKM_rer)[1] <-'Mean_FPKM'

mean(Replicated_Random_FPKM_rer$Mean_FPKM)
mean(RER_genenames$fpkm_mean)

Subsampling_Replicated_RER_FPKM<-ggplot(Replicated_Random_FPKM_rer, aes(x=Mean_FPKM,))+
  geom_histogram(bins = 30,fill='darkgrey')+theme_bw()+
  geom_vline(xintercept = 23.62885, colour='red', linetype='dashed', size=0.6)+
  geom_vline(xintercept = 10.76905, colour='black', linetype='dashed',size=0.6)+
  labs(x='Mean Expression across genes and species (FPKM)', 
       y='count')+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))

Subsampling_Replicated_RER_FPKM

Replicated_Random_FPKM_rer <- as.numeric(Replicated_Random_FPKM_rer$Mean_FPKM)  

wilcox.test(Replicated_Random_FPKM_rer, RER_genenames$fpkm_mean)
wilcox.test(Background$fpkm_mean, RER_genenames$fpkm_mean)

mean(Replicated_Random_FPKM_rer)

#Check that this actually 222 genes.
Unique_Random <- unique(Randomly_sampled_genes_comparison_acceleratedRER$human_alignment_symbol)
length(Unique_Random)

Randomly_sampled_genes_comparison_acceleratedRER<-Randomly_sampled_genes_comparison_acceleratedRER %>%
  mutate(Test='Control')
names(Randomly_sampled_genes_comparison_acceleratedRER)[1] <-'Genenames'

Positive_significant_correlations <-Positive_significant_correlations %>%
  select(Genenames, fpkm_mean) %>%
  mutate(Test='Experimental')


Positive_significant_correlations %>%
  summarise(a=mean(fpkm_mean))


PositiveCorrs_RandomGenes_TestExpression<-do.call("rbind", list(Positive_significant_correlations, 
                        Randomly_sampled_genes_comparison_acceleratedRER))


                        
                        
ggplot(PositiveCorrs_RandomGenes_TestExpression, aes(x=Test, y=fpkm_mean, fill=Test))+
  geom_violin(outlier.shape = NA)+theme_bw()+ylim(0,50)



# Just for show, take the top 30 genes to show
#order by p-value

RER_genenames_plot<-RER_genenames %>%
  slice_min(order_by = -fpkm_mean, n =30)

Negative_plot_RER<-Negative_significant_correlations %>%
  slice_min(order_by = -fpkm_mean, n =50)

Positive_plot_RER<-Positive_significant_correlations %>%
  slice_min(order_by = -fpkm_mean, n =30)

#calculate mean for summary
mean(Background$fpkm_mean)

Negative_RER<-ggplot(Negative_plot_RER, aes(x=fpkm_mean, y=reorder(Genenames, fpkm_mean), colour=P, size=Rho))+
  geom_point(alpha=0.8)+scale_colour_viridis(option='D', discrete = F)+
  geom_vline(xintercept = 10.77793, colour='red', linetype='dashed')+
  theme_bw()+xlab('mean FPKM (Placenta - mammals)')+ylab('Genes')+labs(colour = "P-value")+labs(size = "Rho (Correlation coefficinet)")

Negative_RER

Positive_RER<-ggplot(Positive_plot_RER, aes(x=fpkm_mean, y=reorder(Genenames, fpkm_mean), colour=P, size=Rho))+
  geom_point(alpha=0.8)+scale_colour_viridis(option='D', discrete = F)+
  geom_vline(xintercept = 10.77793, colour='black', linetype='dashed')+
  geom_vline(xintercept = 21.23324, colour='red', linetype='dashed', size=0.6)+
  theme_bw()+xlab('mean FPKM (Placenta - mammals)')+ylab('Genes')+labs(colour = "P-value")+labs(size = "Relative evolutionary rate (Rho)")+
  theme(legend.position = c(0.8, 0.5))+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))

Positive_RER

Top_30_RER<-ggplot(RER_genenames_plot, aes(x=fpkm_mean, y=reorder(Genenames, fpkm_mean), colour=P, size=Rho))+
  geom_point(alpha=0.8)+scale_colour_viridis(option='D', discrete = F)+
  geom_vline(xintercept = 10.77793, colour='black', linetype='dashed')+
  geom_vline(xintercept = 23.62885, colour='red', linetype='dashed', size=0.6)+
  theme_bw()+xlab('mean FPKM (Placenta - mammals)')+ylab('Genes')+labs(colour = "P-value")+labs(size = "Relative evolutionary rate (Rho)")+
  theme(legend.position = c(0.8, 0.5))+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))






library(ggpubr)
tiff("RER_converge_top30_PositiveNegative_Rho_wiAnnots.tiff", units="in", width=10, height=5, res=300)
ggarrange(Positive_RER, Negative_RER, common.legend = T, font.label = list(size= 10, face='bold'), legend = "right" )
dev.off()

tiff("RER_Subsampling_RER_PositiveRER.tiff", units="in", width=11, height=6, res=300)
ggarrange(Subsampling_Replicated_RER_FPKM, Top_30_RER, common.legend = T, font.label = list(size= 10, face='bold'), legend = "right" )
dev.off()


############## ANNOTATIONS ##############
annots=read.gmt("gmtfile.gmt")
stats=getStat(RER_genenames)
annotlist=list(annots)
names(annotlist)="MSigDBpathways"
enrichment=fastwilcoxGMTall(stats, annotlist, outputGeneVals=T, num.g=10)




######### Plotting convergent subs #####

Convergent_subs_null <- read_table2('convergent_genomewide_subs.txt')
Convergent_subs_null_ratio<-Convergent_subs_null %>%
  group_by(Tests) %>%
  mutate(Ratio=Convergent_sites/Sites_Tested)

tiff("Null_convergent_model.tiff", units="in", width=4, height=6, res=300)
ggplot(Convergent_subs_null_ratio, aes(x=Tests, y=Ratio, label=Convergent_sites, fill=Tests))+
  geom_col()+geom_label(size=4)+
  labs(x='Tests conducted',
       y='Convergent sites / Number of sites tested')+
  theme_pubr()
dev.off()



##############################################################################
####################### What proportion of convergent genes also show signals of introgression ##############
##############
#Get RER data with filenames
Viviparity_RER_introgression<-setDT(corViviparity, keep.rownames = "Filenames")
Viviparity_RER_introgression 
#Get filenames for single amino acid convergent changes
Alignments_with_convergent_sites <- read.table('ALIGNMENTS_WITH_CONVERGENT_SITES.txt')
#Now get Introgression data for Goodeids
Goodeid_Introgression_DFOIL <- read.table('Final_Goodeid.dfoil.out.pvals.txt', header=T)
#Now chop off the first bit of the filename for Goodeids
Goodeid_Introgression_DFOIL <- Goodeid_Introgression_DFOIL %>% separate(chrom, c("A", "B", "C", "D"),"/")
Goodeid_Introgression_DFOIL <- Goodeid_Introgression_DFOIL %>% separate(D, c("A", "B"),"\\.")
Goodeid_Introgression_DFOIL <- Goodeid_Introgression_DFOIL %>% separate(A, c("A", "B", "C", "D"),"-")
names(Goodeid_Introgression_DFOIL)[5] <- "scaffold"
names(Goodeid_Introgression_DFOIL)[6] <- "start"
names(Goodeid_Introgression_DFOIL)[7] <- "end"
Goodeid_Introgression_DFOIL$Gene_Identifier<-paste(Goodeid_Introgression_DFOIL$start,Goodeid_Introgression_DFOIL$end)
Goodeid_Introgression_DFOIL<-apply(Goodeid_Introgression_DFOIL,2,function(x)gsub('\\s+', '',x))
# Now prepare the single amino acid changes file
Alignments_with_convergent_sites <- Alignments_with_convergent_sites %>% separate(V1, c("A", "B"),"\\.")
Alignments_with_convergent_sites <- Alignments_with_convergent_sites %>% separate(A, c("A", "B", "C", "D"),"_")
Alignments_with_convergent_sites$D<-sub('.', '', Alignments_with_convergent_sites$D)
Alignments_with_convergent_sites<-as.data.frame(Alignments_with_convergent_sites)
Goodeid_Introgression_DFOIL<-as.data.frame(Goodeid_Introgression_DFOIL)
`%notin%` <- Negate(`%in%`)



#Categorise genes as evidence of convergent evolution and introgression, or not.
Goodeid_Introgression_DFOIL<-Goodeid_Introgression_DFOIL %>%
  mutate(Convergent_amino_acid = case_when(Gene_Identifier %in% Alignments_with_convergent_sites$D ~ "Convergent",
                                       Gene_Identifier %notin% Alignments_with_convergent_sites$D  ~ "Not_Convergent"),
         Introgression = case_when(as.numeric(DFO_Pvalue) < 0.05 | as.numeric(DFI_Pvalue) < 0.05 | as.numeric(DOL_Pvalue) < 0.05 | as.numeric(DIL_Pvalue) < 0.05 ~ "Introgressed",
                                   DFO_Pvalue > 0.05 & DIL_Pvalue > 0.05 & DFI_Pvalue > 0.05 & DOL_Pvalue > 0.05 ~ "Not_introgressed"))

#Convert into 2x2 table
Contigency_Goodeid_Introgression_Convergent<-xtabs(~ Convergent_amino_acid + Introgression, data = Goodeid_Introgression_DFOIL)
#conduct chi-squared test.
chisq.test(Contigency_Goodeid_Introgression_Convergent)
##### Results of chi-squared #######
#X-squared = 8.7598, df = 1, p-value = 0.003079


# Now do the same for Poecilia
Poecilia_Introgression_DFOIL <- read.table('Poecilia.dfoil.out.pvals.txt', header=T)
#Now chop off the first bit of the filename for Poecilia
Poecilia_Introgression_DFOIL <- Poecilia_Introgression_DFOIL %>% separate(chrom, c("A", "B", "C", "D"),"/")
Poecilia_Introgression_DFOIL <- Poecilia_Introgression_DFOIL %>% separate(D, c("A", "B"),"\\.")
Poecilia_Introgression_DFOIL <- Poecilia_Introgression_DFOIL %>% separate(A, c("A", "B", "C", "D"),"-")
names(Poecilia_Introgression_DFOIL)[5] <- "scaffold"
names(Poecilia_Introgression_DFOIL)[6] <- "start"
names(Poecilia_Introgression_DFOIL)[7] <- "end"
Poecilia_Introgression_DFOIL$Gene_Identifier<-paste(Poecilia_Introgression_DFOIL$start,Poecilia_Introgression_DFOIL$end)
Poecilia_Introgression_DFOIL<-apply(Poecilia_Introgression_DFOIL,2,function(x)gsub('\\s+', '',x))

Poecilia_Introgression_DFOIL<-as.data.frame(Poecilia_Introgression_DFOIL)

Poecilia_Introgression_DFOIL<-Poecilia_Introgression_DFOIL %>%
  mutate(Convergent_amino_acid = case_when(Gene_Identifier %in% Alignments_with_convergent_sites$D ~ "Convergent",
                                           Gene_Identifier %notin% Alignments_with_convergent_sites$D  ~ "Not_Convergent"),
         Introgression = case_when(as.numeric(DFO_Pvalue) < 0.05 | as.numeric(DFI_Pvalue) < 0.05 | as.numeric(DOL_Pvalue) < 0.05 | as.numeric(DIL_Pvalue) < 0.05 ~ "Introgressed",
                                   DFO_Pvalue > 0.05 & DIL_Pvalue > 0.05 & DFI_Pvalue > 0.05 & DOL_Pvalue > 0.05 ~ "Not_introgressed"))

#Convert into 2x2 table
Contigency_Poecilia_Introgression_Convergent<-xtabs(~ Convergent_amino_acid + Introgression, data = Poecilia_Introgression_DFOIL)
#conduct chi-squared test.
chisq.test(Contigency_Poecilia_Introgression_Convergent)
###Results for Poecilia
#X-squared = 3.0434, df = 1, p-value = 0.08107




# Now do the same for Poecilia
Xipho_Introgression_DFOIL <- read.table('Xipho.dfoil.out.pvals.txt', header=T)
#Now chop off the first bit of the filename for Poecilia
Xipho_Introgression_DFOIL <- Xipho_Introgression_DFOIL %>% separate(chrom, c("A", "B", "C", "D"),"/")
Xipho_Introgression_DFOIL <- Xipho_Introgression_DFOIL %>% separate(D, c("A", "B"),"\\.")
Xipho_Introgression_DFOIL <- Xipho_Introgression_DFOIL %>% separate(A, c("A", "B", "C", "D"),"-")
names(Xipho_Introgression_DFOIL)[5] <- "scaffold"
names(Xipho_Introgression_DFOIL)[6] <- "start"
names(Xipho_Introgression_DFOIL)[7] <- "end"
Xipho_Introgression_DFOIL$Gene_Identifier<-paste(Xipho_Introgression_DFOIL$start,Xipho_Introgression_DFOIL$end)
Xipho_Introgression_DFOIL<-apply(Xipho_Introgression_DFOIL,2,function(x)gsub('\\s+', '',x))

Xipho_Introgression_DFOIL<-as.data.frame(Xipho_Introgression_DFOIL)

Xipho_Introgression_DFOIL<-Xipho_Introgression_DFOIL %>%
  mutate(Convergent_amino_acid = case_when(Gene_Identifier %in% Alignments_with_convergent_sites$D ~ "Convergent",
                                           Gene_Identifier %notin% Alignments_with_convergent_sites$D  ~ "Not_Convergent"),
         Introgression = case_when(as.numeric(DFO_Pvalue) < 0.05 | as.numeric(DFI_Pvalue) < 0.05 | as.numeric(DOL_Pvalue) < 0.05 | as.numeric(DIL_Pvalue) < 0.05 ~ "Introgressed",
                                   DFO_Pvalue > 0.05 & DIL_Pvalue > 0.05 & DFI_Pvalue > 0.05 & DOL_Pvalue > 0.05 ~ "Not_introgressed"))

#Convert into 2x2 table
Contigency_Xipho_Introgression_Convergent<-xtabs(~ Convergent_amino_acid + Introgression, data = Xipho_Introgression_DFOIL)
#conduct chi-squared test.
chisq.test(Contigency_Xipho_Introgression_Convergent)
#data:  Contigency_Xipho_Introgression_Convergent
#X-squared = 9.3203, df = 1, p-value = 0.002266

Contigency_Xipho_Introgression_Convergent




####What genes are always introgressed and convergent?
Goodeid_Convergence_Introgression_loci<-Goodeid_Introgression_DFOIL %>%
  filter(Introgression=='Introgressed' & Convergent_amino_acid == "Convergent") %>%
  select(scaffold, start, end, Gene_Identifier)


Poecilid_Convergence_Introgression_loci<-Poecilia_Introgression_DFOIL %>%
  filter(Introgression=='Introgressed' & Convergent_amino_acid == "Convergent") %>%
  select(scaffold, start, end, Gene_Identifier)

Xipho_Convergence_Introgression_loci<-Xipho_Introgression_DFOIL %>%
  filter(Introgression=='Introgressed' & Convergent_amino_acid == "Convergent") %>%
  select(scaffold, start, end, Gene_Identifier)


Intersection_Introgressed_loci_Convergent_Loci<-intersect(intersect(Xipho_Convergence_Introgression_loci, 
                    Poecilid_Convergence_Introgression_loci), 
                    Goodeid_Convergence_Introgression_loci)

#Find genes
Intersection_Introgressed_loci_Convergent_Loci$Gene_Identifier_2<-paste(Intersection_Introgressed_loci_Convergent_Loci$scaffold,"-",
                                                                        Intersection_Introgressed_loci_Convergent_Loci$end)
Intersection_Introgressed_loci_Convergent_Loci$transcriptID<-GM_gtf$V9[match(Intersection_Introgressed_loci_Convergent_Loci$Gene_Identifier_2, GM_gtf$Gene_Identifier)]
Intersection_Introgressed_loci_Convergent_Loci$Annotation<-annot$Preferred_name[match(Intersection_Introgressed_loci_Convergent_Loci$transcriptID, annot$query_name)]

intersect(Positive_Selected_Genes_Both_Lineages$Annotation,Intersection_Introgressed_loci_Convergent_Loci$Annotation)




#Quick survey of interesting genes:
convergent_genes_amino_acid_eggnog <- read_tsv('EGGNOG_CONVERGENT_EVOLUTION_FULL_TABLE_ANNOTATIONS.tsv')
mouse_preg_genes<-read_tsv('Genes_mammalian_pregnancy.txt')
immune_genes<-read_tsv('adaptive_immune_genes.txt')
convergent_genes_amino_acid_eggnog$Preferred_name<-str_to_title(convergent_genes_amino_acid_eggnog$Preferred_name)

intersect(mouse_preg_genes$Symbol, convergent_genes_amino_acid_eggnog$Preferred_name)
intersect(immune_genes$Symbol, convergent_genes_amino_acid_eggnog$Preferred_name)


Genes_CNEE[[3]]<-tolower(Genes_CNEE[[3]])
intersect(mouse_pregnant_genes$Symbol, Genes_CNEE$Annotations)



