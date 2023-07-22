library(tidyverse)
library(reshape2)
library(scales)
library(viridis)
library(ggpubr)
library(cowplot)


CNEE_Accelerated_regions <- read.table('All_Accelerated_NonExonic_results_FDRcorrected.fixed.txt', header=T)

CNEE_Accelerated_significant_only<-CNEE_Accelerated_regions %>%
  filter(Accelerated_CNEE_fdr < 0.05)


CNEE_Accelerated_regions_summarised_significant <- CNEE_Accelerated_significant_only %>%
  select(scaffold, start, end, lnlratio, pval,alt_subscale, Accelerated_CNEE_fdr) %>%
  mutate(mid = (start+end)/2)

ggplot(CNEE_Accelerated_regions_summarised_significant, aes(x=scaffold, y=alt_subscale, colour=Accelerated_CNEE_fdr))+
  geom_jitter()+theme_bw()+
  scale_colour_viridis(discrete=F, option='D')+coord_flip()

CNEE_gene_set_Hsap <- read.table("181_GenesNearCNEE_PANTHER_GENE_ONTOLOGY_HUMAN_BACK.txt",
                                       fill = TRUE , header = T)
CNEE_gene_set_Hsap$fold_Enrichment <- as.factor(CNEE_gene_set_Hsap$fold_Enrichment)
CNEE_gene_set_Hsap$GO_biological_process_complete<-gsub("_", " ", CNEE_gene_set_Hsap$GO_biological_process_complete)
top_5_CNEE<-top_n(CNEE_gene_set_Hsap, 5, fold_Enrichment)

glimpse(top_5_CNEE)
top_5_CNEE$fold_Enrichment<-as.numeric(as.character(top_5_CNEE$fold_Enrichment))

tiff("GO_GenesNearCNEE_Panther_181.tiff", units="in", width=7, height=5, res=300)
CNEE_GO<-ggplot(top_5_CNEE , aes(x=fold_Enrichment, y=reorder(GO_biological_process_complete,fold_Enrichment), colour=FDR, size=fold_Enrichment))+
  geom_point()+scale_colour_viridis(option='D', discrete = F)+
  theme_minimal()+xlab('Fold Enrichment')+ylab('')+theme(legend.position = c(0.8, 0.4))+theme(axis.text=element_text(size=12,face='bold'))
CNEE_GO
dev.off()


CNEE_GO


################# Plot altsubscale across genome ##########

CNEE_Accelerated_regions<-CNEE_Accelerated_regions %>%
  mutate(Direction=case_when(
    Accelerated_CNEE_fdr < 0.05 ~ "Significant",
    Accelerated_CNEE_fdr > 0.05 ~ "NonSignificant"))

chrom <- read.csv('All_Chr_Positions_GuppyvsGM.csv', header = F)
CNEE_Accelerated_regions$Chromosome <- chrom$V1[ match(CNEE_Accelerated_regions$scaffold, chrom$V2)]
CNEE_Accelerated_significant_only$Chromosome <- chrom$V1[ match(CNEE_Accelerated_significant_only$scaffold, chrom$V2)]

#sort the chr values because they look crazy
CNEE_Accelerated_regions$Chromosome<-factor(CNEE_Accelerated_regions$Chromosome,levels = c("LG1", 'LG2', 'LG3','LG4','LG5', 'LG6', 'LG7', 'LG8',
                            'LG9', 'LG10', 'LG11', 'LG12','LG13','LG14','LG15', 'LG16', 'LG17',
                            'LG18', 'LG19', 'LG20', 'LG21', 'LG22', 'LG23', 'Scaffold'))


ggplot(CNEE_Accelerated_regions_summarised_significant, aes(x=Chromosome))+
  geom_bar()



#tiff("CNEE_accelerated_across_genome.tiff", units="in", width=10, height=5, res=300)
ggplot(CNEE_Accelerated_regions, aes(x=Chromosome,y=-log10(Accelerated_CNEE_fdr),colour=Direction, alpha=Direction))+
  geom_point()+
  theme_pubr()+
  scale_colour_manual(values = c('darkgrey', 'red'))+
  theme(axis.text.x = element_text(size=0))+
  labs(x='Scaffold', 
       y='Accelerated divergence in foreground branches', 
       colour='FDR')+facet_grid(scaffold ~ ., scales="free")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))

#dev.off()
tiff("CNEE_accelerated_across_genome_test1.tiff", units="in", width=19, height=5, res=300)
ggplot(CNEE_Accelerated_regions, aes(x=Chromosome, y=alt_subscale, colour=Direction, alpha=Direction))+
  geom_jitter(size=1.75,alpha=0.7)+scale_colour_manual(values = c("black", "#FC4E07"))+
  geom_hline(yintercept = 0.0)+
  scale_alpha_manual(guide='none', values = list(Nonsignificant = 0.1, Significant = 1))+
  labs(y=expression('Sequence divergence (alt_subscale)'))+ylim(0,20000)+
  theme_minimal_hgrid()+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
  facet_grid(cols = vars(Chromosome),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Linkage group") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  theme(axis.line.y=element_line(colour = 'grey'))
dev.off()

#####Enrichment of binding motifs 

Binding_motifs <- read.table('AME_MOTIF_BINDING_CNEEs.txt', h=T)
distinct(Binding_motifs, motif_alt_ID)

# Proteins
Protein_domains <- read.csv('Protein_domains_181TFs.csv')

Protein_domains<-top_n(Protein_domains, 10, observed.gene.count)
Protein_domains$observed.gene.count<-as.numeric(Protein_domains$observed.gene.count)
glimpse(Protein_domains)

tiff("Protein_domains_CNEEtfs.tiff", units="in", width=7, height=4, res=300)
ggplot(Protein_domains, aes(x=term.description,y=observed.gene.count))+
  geom_col(fill='red')+theme_pubr()+coord_flip()+
  labs(y='Observed Gene Counts',
       x='Protein domain description')
dev.off()


#Gene ontology for TF binding to CNEEs with significant adjusted p-values.


CNEE_GO_TFs_CNEE_nogenes <- read.csv('GO_Panther_CNEE_putativegenesfilt.csv', header=T)
CNEE_GO_TFs_CNEE_nogenes$Fold.Enrichment<-as.numeric(as.character(CNEE_GO_TFs_CNEE_nogenes$Fold.Enrichment))
#CNEE_GO_TFs_CNEE_nogenes$Gene.ontology<-gsub("<", "", CNEE_GO_TFs_CNEE_nogenes$Gene.ontology)
top_20_CNEE<-top_n(CNEE_GO_TFs_CNEE_nogenes, 20, Fold.Enrichment)
glimpse(top_20_CNEE)

tiff("CNEE_genesremoved_drerio.tiff", units="in", width=7, height=5, res=300)
ggplot(top_20_CNEE, aes(x=Fold.Enrichment, y=reorder(Gene.ontology,Fold.Enrichment), colour=adj.Pvalue, size=Fold.Enrichment))+
  geom_point()+scale_colour_viridis(option='D', discrete = F)+
  theme_minimal()+xlab('Fold Enrichment')+ylab('Biological process - Zebrafish background')+
  labs(colour='adjusted p-value',
       size= 'Fold Enrichment')
dev.off()


####################################################################
############### Number of CNEEs near genes ########################

Genes_CNEE <- read.table('Genes_NumberofAccCNEESperGene.txt', header=F)
Genes_CNEE$Annotation<-annot$Preferred_name[match(Genes_CNEE$V2, annot$query_name)]
Genes_CNEE$Description<-annot$X.4[match(Genes_CNEE$V2, annot$query_name)]
colnames(Genes_CNEE) <- c('Number_of_CNEE', 'transcript_ID', 'Annotations')

Genes_More_Than_One_CNEE <- Genes_CNEE %>%
  filter(Number_of_CNEE > 1)

tiff("CNEE_NumberofCNEENextToGene.tiff", units="in", width=6, height=4, res=300)
ggplot(Genes_CNEE, aes(x=Number_of_CNEE))+ylim(0,200)+
  geom_bar()+xlab('Number of accelerated CNE near gene')+
  geom_text(stat='count', aes(label=..count..), hjust=-0.3)+scale_fill_manual('black')+
  theme_pubr()+coord_flip()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))
dev.off()




########## Test mammalian expression of non-coding genes
Genes_CNEE_uppercase<-mutate_all(Genes_CNEE, .funs = toupper)
Genes_CNEE_uppercase$fpkm_mean <- Background$fpkm_mean[ match(Genes_CNEE_uppercase$Annotations, Background$human_alignment_symbol)]
