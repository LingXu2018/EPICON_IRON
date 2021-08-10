setwd("~/Downloads/Iron_MS/")
library("phyloseq")
library("ggplot2")
library("scales")
library("grid")
theme_set(theme_bw())
library("DESeq2")
library("ape")
library("vegan")
library("data.table")
library("RColorBrewer")
library(colorRamps)
library("svglite")
library(VennDiagram)

biom_file <- "mag_input_biom.biom"
map_file <- "metadata_MAG.txt"  
set.seed(1)
biomot=import_biom(biom_file,parseFunction = parse_taxonomy_greengenes)
bmsd = import_qiime_sample_data(map_file)
Bushman = merge_phyloseq(biomot,bmsd)
##########################################
############figure 1e-1f, supp.4e-4f###########
#########################################
options(scipen=200)
b1_phylum <- tax_glom(Bushman, taxrank="Phylum")

rt430=subset_samples(b1_phylum,Genotype=="RT430")
bt642=subset_samples(b1_phylum,Genotype=="T642")

c = subset_samples(b1_phylum,Treatment=="Control")
x <- c("TP3","TP4","TP8","TP9","TP10","TP11")
rm <- {}
control <- {}
for(i in x){
  TP = subset_samples(c, Timepoint == i)
  merged_r = merge_samples(TP, "SampleType")
  merged_r_m = transform_sample_counts(merged_r , function(x) 100 * x/sum(x))
  physeqdf <- psmelt(merged_r_m)
  temp <- c("Soil","Rhizosphere","Root")
  physeqdf$Sample <- factor(physeqdf$Sample, levels = c(temp))
  n = nrow(physeqdf[order(physeqdf$Phylum),])
  x = sample(i,n,replace=TRUE)
  ST <- sample("Control",n,replace=TRUE)
  physeqdf[order(physeqdf$Phylum),]$Timepoint <- x 
  physeqdf[order(physeqdf$Phylum),]$Treatment <- ST
  physeqdf$Phylum <- as.character(physeqdf$Phylum)
  #physeqdf$Phylum[30:129] <- "Other"
  control <- rbind(control,physeqdf[order(physeqdf$Phylum),])
}

control <- data.frame(control)


pre2 <- {}
pre = subset_samples(bt642,Treatment=="Pre_flowering")
y <- c("TP3","TP4","TP8","TP9","TP10")
for(i in y){
  TP = subset_samples(pre, Timepoint == i)
  merged_r = merge_samples(TP, "SampleType")
  merged_r_m = transform_sample_counts(merged_r , function(x) 100 * x/sum(x))
  physeqdf <- psmelt(merged_r_m)
  temp <- c("Soil","Rhizosphere","Root")
  physeqdf$Sample <- factor(physeqdf$Sample, levels = c(temp))
  n = nrow(physeqdf[order(physeqdf$Phylum),])
  x = sample(i,n,replace=TRUE)
  ST <- sample("Pre_flowering",n,replace=TRUE)
  physeqdf[order(physeqdf$Phylum),]$Timepoint <- x 
  physeqdf[order(physeqdf$Phylum),]$Treatment <- ST
  physeqdf$Phylum <- as.character(physeqdf$Phylum)
  #physeqdf$Phylum[30:129] <- "Other"
  pre2 <- rbind(pre2,physeqdf[order(physeqdf$Phylum),])
}

add2 <- subset(control,Treatment=="Control")
add3 <- subset(add2,Timepoint=="TP3"|Timepoint=="TP4"|Timepoint=="TP8"|Timepoint=="TP9")
add3$Treatment[add3$Treatment == "Control"] <- "Post_flowering"

p <- {}
y <- c("TP10","TP11")
post = subset_samples(bt642,Treatment=="Post_flowering")
for(i in y){
  TP = subset_samples(post, Timepoint == i)
  merged_r = merge_samples(TP, "SampleType")
  merged_r_m = transform_sample_counts(merged_r , function(x) 100 * x/sum(x))
  physeqdf <- psmelt(merged_r_m)
  temp <- c("Soil","Rhizosphere","Root")
  physeqdf$Sample <- factor(physeqdf$Sample, levels = c(temp))
  n = nrow(physeqdf[order(physeqdf$Phylum),])
  x = sample(i,n,replace=TRUE)
  ST <- sample("Post_flowering",n,replace=TRUE)
  physeqdf[order(physeqdf$Phylum),]$Timepoint <- x 
  physeqdf[order(physeqdf$Phylum),]$Treatment <- ST
  physeqdf$Phylum <- as.character(physeqdf$Phylum)
  p <- rbind(p,physeqdf[order(physeqdf$Phylum),])
}
rm <- rbind(control,pre2,add3,p)

sum <- tapply(rm$Abundance,rm$Phylum, FUN=sum)
sum <- data.frame(sum)
sum <- sort(sum$sum,decreasing = T)
list <- c("Acidobacteria","Actinobacteria","Bacteroidetes","Chloroflexi","Cyanobacteria",
          "Firmicutes","Gemmatimonadetes","Nitrospirae",
          "Planctomycetes","Proteobacteria","Tenericutes","Saccharibacteria","Verrucomicrobia","Unclassified")
a <- rm[rm$Phylum %in% list,]
b <- rm[!(rm$Phylum %in% list),]
b$Phylum <- "Other"
rm <- rbind(a,b)
rm$Treatment<-factor(rm$Treatment,levels=c("Control","Pre_flowering","Post_flowering"))
rm$Timepoint<-factor(rm$Timepoint, levels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17"))
rm$Phylum[rm$Phylum=="Saccharibacteria"] <- "Other"

col_SampleByTreatment<-c("Acidobacteria"="#4d5b6a","Actinobacteria"="#6a90b4","
                         Armatimonadetes"="#599ec4","Bacteroidetes"="#a1c5c8",
                         "Chloroflexi"="#84b59c","Cyanobacteria"="#427e83",
                         "Firmicutes"="#7997a1","Gemmatimonadetes"="#919f70",
                         "Nitrospirae"="#686a47","Planctomycetes"="#c8bd83",
                         "Proteobacteria"="#cb9e59","SPAM"="#ecc363","Tenericutes"="#c57b67","TM7"="#d65f3d",
                         "Verrucomicrobia"="#a04344","Other"="#bf7a8f","Saccharibacteria"="#bf7a8f")




rm$Phylum <- factor(rm$Phylum,levels=c("Actinobacteria","Chloroflexi","Bacteroidetes","Proteobacteria",
                                       "Acidobacteria","Gemmatimonadetes",
                                       "Nitrospirae","Verrucomicrobia","Other"))

plot <- ggplot(rm, aes(x=Timepoint, y=Abundance, fill=Phylum)) + 
  geom_bar(stat = "identity", color = "NA")+facet_grid(Sample~Treatment,scales = "free_x",space="free_x")+ scale_fill_manual(values = col_SampleByTreatment) +
  theme(axis.text.x=element_text(size=11,color="black",angle=90), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=14,face="bold"),text=element_text(size=14,face="bold"))+
  ylab("Relative Abundance")+xlab("")

ggsave(filename = "RA-bins-newcolor.jpg",plot = plot,width=10,height=6)


###########################supplementary figure 7###############
metadata_bushman=data.frame(sample_data(Bushman))
metadata_bushman$Location<-"Kearney"
sample_data(Bushman)=metadata_bushman

data=merge_samples(Bushman,"Location")
relative_abundance = transform_sample_counts(data , function(x) 100 * x/sum(x))
bin=data.frame(t(otu_table(relative_abundance)))
tax=data.frame(tax_table(relative_abundance))
rm=cbind(bin,tax)

col_SampleByTreatment<-c("Acidobacteria"="#4d5b6a","Actinobacteria"="#6a90b4","
                         Armatimonadetes"="#599ec4","Bacteroidetes"="#a1c5c8",
                         "Chloroflexi"="#84b59c","Cyanobacteria"="#427e83",
                         "Firmicutes"="#7997a1","Gemmatimonadetes"="#919f70",
                         "Nitrospirae"="#686a47","Planctomycetes"="#c8bd83",
                         "Proteobacteria"="#cb9e59","SPAM"="#ecc363",
                         "Tenericutes"="#c57b67","TM7"="#d65f3d",
                         "Verrucomicrobia"="#a04344","Other"="#bf7a8f",
                         "Saccharibacteria"="#bf7a8f")


rm$Phylum <- factor(rm$Phylum,levels=c("Actinobacteria","Chloroflexi","Bacteroidetes","Proteobacteria",
                                       "Acidobacteria","Gemmatimonadetes",
                                       "Nitrospirae","Verrucomicrobia","Saccharibacteria"))
#supplementary figure 7b
##########Relative abundance for each bin in 52 samples
p1=ggplot(rm, aes(x=Phylum, y=Kearney, fill=Phylum)) + 
  geom_bar(stat = "identity", color = "black")+#facet_grid(Sample~Treatment,scales = "free_x",space="free_x")+ 
  scale_fill_manual(values = col_SampleByTreatment) +
  theme(axis.text.x=element_text(size=11,color="black",angle=90), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=14),text=element_text(size=14))+xlab("")+
  ylab("Relative Abundance")
ggplot(ggsave(filename = "relative_abundance_bins_52.pdf",plot = p1,width=5,height=5))

#supplementary figure 7a
###########Bin numbers for each phylum
p2=ggplot(data=rm, aes(x=Phylum,fill=Phylum)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-1)+
  scale_fill_manual(values = col_SampleByTreatment)+
  theme(axis.text.x=element_text(size=11,color="black",angle=90), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=14),text=element_text(size=14))+
  ylab("Bin Count")+xlab("")
ggplot(ggsave(filename = "Bin_count_52.pdf",plot = p2,width=5,height=5))

#####figure2a-2b, supp.figure8,shotgun function COGs
library(data.table)
library(ggplot2)
library(RColorBrewer)
setwd("/Users/lingxu/Downloads/Shotgun_bins_materials/Two_groups_of_genome_comparison/")
cog_detail=read.table("cog_detail.txt")
colnames(cog_detail)=c("COG_ID","Characters")

funct <- read.table("function_characters.txt")
colnames(funct)<-c("Characters","Annotation")

#Supplementary Fig. 8
#########################################
#Actinobacteria vs Actinobacteria
#########################################
##########actino_vs_actino_enriched_D
aae=read.csv("actino_vs_actino_enriched_drought.txt",sep = "\t",skip = 0, header = FALSE,
             comment.char = "", check.names = FALSE)
colnames(aae)=c("COG_ID","Description")
mer_aae<-merge(aae,cog_detail,by="COG_ID")
mer_aae=merge(mer_aae,funct,by="Characters")

mer_aae_new=as.data.frame(table(mer_aae$Characters))
colnames(mer_aae_new)<-c("Characters","AAE")

all=as.data.frame(table(cog_detail$Characters))
colnames(all)=c("Characters","all_cog")
combine_AAE=merge(mer_aae_new,all,by="Characters")
combine_AAE$hyper=phyper(combine_AAE$AAE-1,combine_AAE$all_cog,4384-combine_AAE$all_cog,235,lower.tail=FALSE)
combine_AAE$aae_ratio=combine_AAE$AAE/235
combine_AAE$all_ratio=combine_AAE$all_cog/4384
combine_AAE$final_ratio=combine_AAE$aae_ratio/combine_AAE$all_ratio
characters_aae=merge(funct,combine_AAE,by="Characters")
characters_aae$Status="enriched_AA"
characters_aae=characters_aae[,c(1,2,5,8,9)]
colnames(characters_aae)=c("Characters","Annotation","hyper","Final_ratio","Status")

##########actino_vs_actino_depletion_drought
aad=read.csv("actino_vs_actino_depletion_drought.txt",sep = "\t",skip = 0, header = FALSE,
             comment.char = "", check.names = FALSE)
colnames(aad)=c("COG_ID","Description")
mer_aad<-merge(aad,cog_detail,by="COG_ID")
mer_aad=merge(mer_aad,funct,by="Characters")
#write.csv(mer_aad,"actino_vs_actino_drought_depleted_COG.csv")

mer_aad_new=as.data.frame(table(mer_aad$Characters))
colnames(mer_aad_new)<-c("Characters","AAD")

all=as.data.frame(table(cog_detail$Characters))
colnames(all)=c("Characters","all_cog")
combine_AAD=merge(mer_aad_new,all,by="Characters")
combine_AAD$hyper=phyper(combine_AAD$AAD-1,combine_AAD$all_cog,4384-combine_AAD$all_cog,213,lower.tail=FALSE)
combine_AAD$aad_ratio=combine_AAD$AAD/213
combine_AAD$all_ratio=combine_AAD$all_cog/4384
combine_AAD$final_ratio=combine_AAD$aad_ratio/combine_AAD$all_ratio
characters_AAD=merge(funct,combine_AAD,by="Characters")
characters_AAD$Status="depleted_AA"
#write.csv(characters_AAD,"actino_vs_actino_drought_depleted_high_level_COG.csv")
characters_AAD=characters_AAD[,c(1,2,5,8,9)]
colnames(characters_AAD)=c("Characters","Annotation","hyper","Final_ratio","Status")

final=rbind(characters_AAD,characters_aae)

final$Annotation <- factor(final$Annotation, levels=c("Carbohydrate_transport_and_metabolism",
                                                      "Transcription",
                                                      "Inorganic_ion_transport_and_metabolism",
                                                      "Amino_acid_transport_and_metabolism",
                                                      "Coenzyme_transport_and_metabolism",
                                                      "Replication_recombination_and_repair",
                                                      "Energy_production_and_conversion",
                                                      "Posttranslational_modification_protein_turnover_chaperones",
                                                      "RNA_processing_and_modification",
                                                      "Lipid_transport_and_metabolism",
                                                      "Signal_transduction_mechanisms",
                                                      "Secondary_metabolites_biosynthesis_transport_and_catabolism",
                                                      "Cell_cycle_control_cell_division_chromosome_partitioning",
                                                      "Defense_mechanisms",
                                                      "General_function_prediction_only",
                                                      "Nucleotide_transport_and_metabolism",
                                                      "Cell_wall/membrane/envelope_biogenesis",
                                                      "Intracellular_trafficking_secretion_and_vesicular_transport",
                                                      "Mobilome_prophages_transposons",
                                                      "Translation_ribosomal_structure_and_biogenesis",
                                                      "Function_unknown",
                                                      "Chromatin_structure_and_dynamics",
                                                      "Cell_motility",
                                                      "Extracellular_structures",
                                                      "Cytoskeleton"))


final_new<-subset(final, Annotation!="Cytoskeleton" & Annotation!="Chromatin_structure_and_dynamics"&
                    Annotation!="RNA_processing_and_modification")

p=ggplot(final_new, aes(x=Final_ratio, y=Annotation,color=Status)) +facet_wrap(~Status,scales="free_x")+xlab("Ratio")+
  scale_colour_manual(name="",values = c("#004445","#6FB98F","#9B4F0F","#C99E10"))+
  scale_shape_manual(name="",values=c(1,1,19,19))+geom_point(size = 4,stroke = 1)+theme_bw()+
  ylab("")+
  geom_vline(aes(xintercept=1.5),colour="Black",size=1,linetype="dashed")+xlab("Enrichment ratio")+
  theme(axis.text.x=element_text(size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=10))+
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(strip.text = element_text(size=11,face="bold"))+ggtitle("Actino_vs_Actino")

ggsave(filename = "High_level_function_AA_COG.pdf",plot = p,width=12,height=6)


##########actino_vs_actino_depletion_drought
aad=read.csv("actino_vs_proteo_depleted_drought.txt",sep = "\t",skip = 0, header = FALSE,
             comment.char = "", check.names = FALSE)
colnames(aad)=c("COG_ID","Description")
mer_aad<-merge(aad,cog_detail,by="COG_ID")
mer_aad=merge(mer_aad,funct,by="Characters")
#write.csv(mer_aad,"actino_vs_actino_drought_depleted_COG.csv")

mer_aad_new=as.data.frame(table(mer_aad$Characters))
colnames(mer_aad_new)<-c("Characters","AAD")

all=as.data.frame(table(cog_detail$Characters))
colnames(all)=c("Characters","all_cog")
combine_AAD=merge(mer_aad_new,all,by="Characters")
combine_AAD$hyper=phyper(combine_AAD$AAD-1,combine_AAD$all_cog,4384-combine_AAD$all_cog,213,lower.tail=FALSE)
combine_AAD$aad_ratio=combine_AAD$AAD/213
combine_AAD$all_ratio=combine_AAD$all_cog/4384
combine_AAD$final_ratio=combine_AAD$aad_ratio/combine_AAD$all_ratio
characters_AAD=merge(funct,combine_AAD,by="Characters")
characters_AAD$Status="depleted_AA"
#write.csv(characters_AAD,"actino_vs_actino_drought_depleted_high_level_COG.csv")
characters_AAD=characters_AAD[,c(1,2,5,8,9)]
colnames(characters_AAD)=c("Characters","Annotation","hyper","Final_ratio","Status")

final=rbind(characters_AAD,characters_aae)

final$Annotation <- factor(final$Annotation, levels=c("Carbohydrate_transport_and_metabolism",
                                                      "Transcription",
                                                      "Inorganic_ion_transport_and_metabolism",
                                                      "Amino_acid_transport_and_metabolism",
                                                      "Coenzyme_transport_and_metabolism",
                                                      "Replication_recombination_and_repair",
                                                      "Energy_production_and_conversion",
                                                      "Posttranslational_modification_protein_turnover_chaperones",
                                                      "RNA_processing_and_modification",
                                                      "Lipid_transport_and_metabolism",
                                                      "Signal_transduction_mechanisms",
                                                      "Secondary_metabolites_biosynthesis_transport_and_catabolism",
                                                      "Cell_cycle_control_cell_division_chromosome_partitioning",
                                                      "Defense_mechanisms",
                                                      "General_function_prediction_only",
                                                      "Nucleotide_transport_and_metabolism",
                                                      "Cell_wall/membrane/envelope_biogenesis",
                                                      "Intracellular_trafficking_secretion_and_vesicular_transport",
                                                      "Mobilome_prophages_transposons",
                                                      "Translation_ribosomal_structure_and_biogenesis",
                                                      "Function_unknown",
                                                      "Chromatin_structure_and_dynamics",
                                                      "Cell_motility",
                                                      "Extracellular_structures",
                                                      "Cytoskeleton"))

p=ggplot(final, aes(x=Final_ratio, y=Annotation,color=Status)) +facet_wrap(~Status)+xlab("Ratio")+
  scale_colour_manual(name="",values = c("#004445","#6FB98F","#9B4F0F","#C99E10"))+
  scale_shape_manual(name="",values=c(1,1,19,19))+geom_point(size = 4,stroke = 1)+theme_bw()+
  ylab("")+
  geom_vline(aes(xintercept=1.5),colour="Black",size=1,linetype="dashed")+xlab("Enrichment ratio")+
  theme(axis.text.x=element_text(size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=10))+
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(strip.text = element_text(size=11,face="bold"))+ggtitle("Actino_vs_Proteo")

ggsave(filename = "High_level_function_AP_COG.pdf",plot = p,width=12,height=8)
