#Expression data for placental mammals - Look at and plot expression data for candidate genes.
library(tidyverse)
library(reshape2)
library(scales)
library(viridis)
#First import expression data
A_fusciceps<-read.csv("CSV_placenta_mammals/Ateles_fusciceps.csv", header = T, fill=T, row.names=NULL)
H_sapiens<-read.csv("CSV_placenta_mammals/Homo_sapiens.csv", header = T, fill=T, row.names=NULL)
B_taurus<-read.csv("CSV_placenta_mammals/Bos_taurus.csv", header = T, fill=T, row.names=NULL)
C_familiaris<-read.csv("CSV_placenta_mammals/Canis_familiaris.csv", header = T, fill=T, row.names=NULL)
D_novemcinctus<-read.csv("CSV_placenta_mammals/Dasypus_novemcinctus.csv", header = T, fill=T, row.names=NULL)
E_callabus<-read.csv("CSV_placenta_mammals/Equus_caballus.csv", header = T, fill=T, row.names=NULL)
L_africana<-read.csv("CSV_placenta_mammals/Loxodonta_africana.csv", header = T, fill=T, row.names=NULL)
M_domestica<-read.csv("CSV_placenta_mammals/Monodelphis_domestica.csv", header = T, fill=T, row.names=NULL)
M_musculus<-read.csv("CSV_placenta_mammals/Mus_musculus.csv", header = T, fill=T, row.names=NULL)
O_aries<-read.csv("CSV_placenta_mammals/Ovis_aries.csv", header = T, fill=T, row.names=NULL)
P_paniscus<-read.csv("CSV_placenta_mammals/Pan_paniscus.csv", header = T, fill=T, row.names=NULL)
S_carmeli<-read.csv("CSV_placenta_mammals/Spalax_carmeli.csv", header = T, fill=T, row.names=NULL)
S_galili<-read.csv("CSV_placenta_mammals/Spalax_galili.csv", header = T, fill=T, row.names=NULL)
S_scrofa<-read.csv("CSV_placenta_mammals/Sus_scrofa.csv", header = T, fill=T, row.names=NULL)

#Next thing to do is simplify the datasets by only keeping rows that contain gene symbol and expression data
#Humans
H_sapiens<-select(H_sapiens, -c(row.names, X))
H_sapiens <- H_sapiens[!apply(is.na(H_sapiens) | H_sapiens == "", 1, all),]
#Spider-headed monkey
A_fusciceps<-select(A_fusciceps, -c(row.names, X))
A_fusciceps <- A_fusciceps[!apply(is.na(A_fusciceps) | A_fusciceps == "", 1, all),]
#Cow
B_taurus<-select(B_taurus, -c(row.names, X))
B_taurus <- B_taurus[!apply(is.na(B_taurus) | B_taurus == "", 1, all),]
#Dog
C_familiaris<-select(C_familiaris, -c(row.names, X))
C_familiaris <- C_familiaris[!apply(is.na(C_familiaris) | C_familiaris == "", 1, all),]
#Armadillo - remove brackets from this.
D_novemcinctus<-select(D_novemcinctus, -c(Gene.symbol, X))
D_novemcinctus <- D_novemcinctus[!apply(is.na(D_novemcinctus) | D_novemcinctus == "", 1, all),]
D_novemcinctus<- D_novemcinctus %>% mutate(Human.Symbol = str_replace_all(Human.Symbol, "\\*|\\(|\\)", ""))
#Horse
E_callabus<-select(E_callabus, -c(row.names, X))
E_callabus<- E_callabus[!apply(is.na(E_callabus) | E_callabus == "", 1, all),]
#Elephant - need to figure this one out.
L_africana<-select(L_africana, -c(Gene.symbol, X, Human.Name))
L_africana <- L_africana[!apply(is.na(L_africana) | L_africana == "", 1, all),]
L_africana<- L_africana %>% mutate(Human.Name = str_replace_all(Human.Name, "\\*|\\(|\\)", ""))
#Opposum


#Mouse-messed up.
M_musculus<-select(M_musculus, -c(row.names, X))
M_musculus<- M_musculus[!apply(is.na(M_musculus) | M_musculus == "", 1, all),]

#Sheep
O_aries<-select(O_aries, -c(row.names, X, Human.Symbol))
O_aries<- O_aries[!apply(is.na(O_aries) | O_aries == "", 1, all),]
View(O_aries)
#Bonobo 
P_paniscus <-select(P_paniscus, -c(row.names, X, Human.Symbol))
P_paniscus <- P_paniscus[!apply(is.na(P_paniscus) | P_paniscus== "", 1, all),]

#S_carmeli
S_carmeli <-select(S_carmeli, -c(row.names, X, Human.Symbol))
S_carmeli <- S_carmeli[!apply(is.na(S_carmeli) | S_carmeli== "", 1, all),]
S_carmeli<- S_carmeli %>% mutate(Human.Symbol = str_replace_all(Human.Symbol, "\\*|\\(|\\)", ""))

#S_galili
S_galili <-select(S_galili, -c(X.1, X, Human.Symbol))
S_galili <- S_galili[!apply(is.na(S_galili) | S_galili== "", 1, all),]
S_galili<- S_galili %>% mutate(Human.Symbol = str_replace_all(Human.Symbol, "\\*|\\(|\\)", ""))

#S_carmeli
S_scrofa <-select(S_scrofa, -c(row.names, X, Human.Symbol))
S_scrofa <- S_carmeli[!apply(is.na(S_scrofa) | S_scrofa== "", 1, all),]
S_scrofa<- S_scrofa %>% mutate(Human.Symbol = str_replace_all(Human.Symbol, "\\*|\\(|\\)", ""))

#Ok - now extract the genes that were found in our convergence analysis.
convergent_genes<-read.table('GENENAMES_59_GENES_CODONS.txt', header=F)
RER_genes<-read.table('RER_rapid_proteins_genenames.txt', header=F)



#Ok I manually checked and compiled list of genes and there expression across mammalian placenta.
#Read this table in.
Placenta_mammalian_expression<-read.csv("Placenta_expression_convergent_genes.csv")
#Melt the csv so it's more useable
Melt_Placenta_mam_exp<-melt(Placenta_mammalian_expression, id.vars = 'Genes')
#Now plot! 
sapply(Melt_Placenta_mam_exp, class)
Melt_Placenta_mam_exp$value<-as.numeric(as.character(Melt_Placenta_mam_exp$value))
Melt_Placenta_mam_exp$variable<-gsub( "_", ".", as.character(Melt_Placenta_mam_exp$variable))

png('Placental_expression_mammals_59_genes.png', width=10, height=8, units='in', res=320)
ggplot(Melt_Placenta_mam_exp, aes(x=Genes,y=value,colour=variable))+
  geom_point(size=2.6, alpha=0.6)+xlab('Genes')+ylab('FPKM (Placenta)')+coord_flip()+
  theme_bw()+labs(colour="Species")
dev.off()

#read in mean transcript data per species from mammalian core transcriptome paper
Mean_core_transcript <- read.table('all_species_mean_fpkm.txt')
#match genes from 59 geneset and 32
Mean_core_convergent<-Mean_core_transcript[(Mean_core_transcript$human_alignment_symbol%in%convergent_genes$V1),]
Mean_RER<-Mean_core_transcript[(Mean_core_transcript$human_alignment_symbol%in%RER_genes$V1),]



#remove weird duplicates from 59 genes
nodups_Mean_core_convergent<-Mean_core_convergent %>% distinct(human_alignment_symbol, mean_fpkm, species, human_name ,.keep_all = TRUE)
nodups_Mean_core_convergent<-nodups_Mean_core_convergent[!is.na(nodups_Mean_core_convergent$human_name), ]
nodups_Mean_core_convergent<-filter(nodups_Mean_core_convergent, human_alignment_symbol != "CASR")

RER_nodups<-Mean_RER %>% distinct(human_alignment_symbol, mean_fpkm, species, human_name ,.keep_all = TRUE)


###### SUMMARISE mean_fpkm by species and gene 
FPKM_grouped_by_gene_and_species<-nodups_Mean_core_convergent %>%
  group_by(species, human_alignment_symbol) %>%
  summarise(mean_fpkm)



####THIS IS THE ONE!
tiff("59_genes_placenta_expression_mamals.tiff", units="in", width=11, height=6.5, res=300)
ggplot(FPKM_grouped_by_gene_and_species, aes(x=human_alignment_symbol,y=mean_fpkm))+
  xlab('Genes')+ylab('FPKM (Placenta)')+coord_flip()+
  geom_point(aes(colour=species))+geom_boxplot(outlier.shape = NA)+
  theme_bw()+labs(colour="Species")+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                                  labels = trans_format("log10", math_format(10^.x)))+
  scale_colour_viridis(discrete=T, option='D')+geom_hline(yintercept = 1, colour='red', linetype='dashed')
dev.off()


#########-----------
mammal_species_aggregate<-aggregate(nodups_Mean_core_convergent[, 14], list(nodups_Mean_core_convergent$human_name), mean)
options(scipen = 999)

mammal_species_aggregate<-filter(mammal_species_aggregate, x > 1)


###### PLOT RER genes but first subset

Background<-Mean_core_transcript %>% 
  group_by(human_alignment_symbol) %>%
  summarise(fpkm_mean=mean(FPKM))

Summary_rer_nodups<-RER_nodups %>% 
  group_by(human_alignment_symbol) %>%
  summarise(fpkm_mean=mean(FPKM))

Summary_rer_larger_than_1<-Summary_rer_nodups %>%
  filter(fpkm_mean >1)



  

######## POECILIOPSIS EXPRESSION #########
 
Poeciliopsis <- read.csv('Poeciliopsis_expression_limited_gene_expression.csv')
Both_pret_pturr<-Poeciliopsis[(Poeciliopsis$both%in%convergent_genes$V1),]


#what's the overlap between mammals and poeciliopsis with non-negigible expression
Mammals_poeciliopsis_both<-mammal_species_aggregate[(mammal_species_aggregate$Group.1%in%Both_pret_pturr$both),]
Mammals_poeciliopsis_tur<-mammal_species_aggregate[(mammal_species_aggregate$Group.1%in%Pturr$P.tur_only),]
Mammals_poeciliopsis_ret<-mammal_species_aggregate[(mammal_species_aggregate$Group.1%in%Pret$P.ret_only),]


######### 
#RER CONVERGENCE VS EXPRESSION
#Any difference between genes with higher expression compared to genes with lower expression.
RER_all_genes <- read.csv('RER_CORRELATIONS_WITH_GENENAMES.csv')
