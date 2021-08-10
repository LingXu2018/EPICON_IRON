setwd("/Users/lingxu/Downloads/Iron_MS/")
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

Bushman <- readRDS("Bushman_shotgun_whole.rds")

sample_data(Bushman)$SampleType<-factor(sample_data(Bushman)$SampleType, levels=c("Soil","Rhizosphere"))
sample_data(Bushman)$Treatment<-factor(sample_data(Bushman)$Treatment, levels=c("Control","Pre_flowering","Post_flowering"))
sample_data(Bushman)$Timepoint<-factor(sample_data(Bushman)$Timepoint, levels=c("TP3","TP4","TP8","TP9","TP10","TP11"))

Bushman = subset_taxa(Bushman, Phylum != "Streptophyta")

############figure1c-1d, supp.4c-4d
####################Combined RT430 & BT642 composition profile
Phylum_rhi <- readRDS("shotgun_whole_rhi_Phylum_Bacteria.rds")

x <- c("TP3","TP4","TP8","TP9","TP10","TP11")
rm <- {}
for(i in x){
  TP = subset_samples(Phylum_rhi, Timepoint == i)
  merged_r = merge_samples(TP, "Treatment")
  merged_r_m = transform_sample_counts(merged_r , function(x) 100 * x/sum(x))
  physeqdf <- psmelt(merged_r_m)
  #temp <- c("BT642Control","RT430Control","BT642Pre_flowering","RT430Pre_flowering","BT642Post_flowering","RT430Post_flowering")
  temp <- c("Control","Pre_flowering","Post_flowering")
  physeqdf$Sample <- factor(physeqdf$Sample, levels = c(temp))
  n = nrow(physeqdf[order(physeqdf$Phylum),])
  x = sample(i,n,replace=TRUE)
  ST <- sample("Rhizosphere",n,replace=TRUE)
  physeqdf[order(physeqdf$Phylum),]$Timepoint <- x 
  physeqdf[order(physeqdf$Phylum),]$SampleType <- ST
  rm <- rbind(rm,physeqdf[order(physeqdf$Phylum),])
}

Phylum_soil <- readRDS("shotgun_whole_soil_Phylum_Bacteria.rds")

x <- c("TP3","TP4","TP8","TP9","TP10","TP11")
rm2 <- {}
for(i in x){
  TP = subset_samples(Phylum_soil, Timepoint == i)
  merged_r = merge_samples(TP, "Treatment")
  merged_r_m = transform_sample_counts(merged_r , function(x) 100 * x/sum(x))
  physeqdf <- psmelt(merged_r_m)
  temp <- c("Control","Pre_flowering","Post_flowering")
  physeqdf$Sample <- factor(physeqdf$Sample, levels = c(temp))
  n = nrow(physeqdf[order(physeqdf$Phylum),])
  x = sample(i,n,replace=TRUE)
  ST <- sample("Soil",n,replace=TRUE)
  physeqdf[order(physeqdf$Phylum),]$Timepoint <- x 
  physeqdf[order(physeqdf$Phylum),]$SampleType <- ST
  rm2 <- rbind(rm2,physeqdf[order(physeqdf$Phylum),])
}

rm3 <- rbind(rm,rm2)

sum <- tapply(rm3$Abundance,rm3$Phylum, FUN=sum)
sum <- data.frame(sum)
sum <- sort(sum$sum,decreasing = T)
#list <- rownames(sum[1:10])
list <- c("Acidobacteria","Actinobacteria","Bacteroidetes","Chloroflexi","Cyanobacteria",
          "Firmicutes","Gemmatimonadetes","Nitrospirae",
          "Planctomycetes","Proteobacteria","Tenericutes","TM7","Verrucomicrobia")

a <- rm3[rm3$Phylum %in% list,]
b <- rm3[!(rm3$Phylum %in% list),]
b$Phylum <- "Other"
rm <- rbind(a,b)
rm$Treatment<-factor(rm$Treatment,levels=c("Control","Pre_flowering","Post_flowering"))
rm$Timepoint<-factor(rm$Timepoint, levels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17"))

# col_SampleByTreatment<-c("Acidobacteria"="thistle","Actinobacteria"="lightsalmon2","Armatimonadetes"="mediumpurple1","Bacteroidetes"="cornsilk","Chloroflexi"="lightskyblue","Cyanobacteria"="lightpink","Firmicutes"="lightblue1","Gemmatimonadetes"="darkseagreen","Nitrospirae"="mistyrose2",
#                          "Planctomycetes"="gray69","Proteobacteria"="cadetblue3","SPAM"="darkorchid1","Tenericutes"="darkolivegreen2","TM7"="yellow","Verrucomicrobia"="slategray1","Other"="plum")

col_SampleByTreatment<-c("Acidobacteria"="#4d5b6a","Actinobacteria"="#6a90b4","Armatimonadetes"="#599ec4","Bacteroidetes"="#a1c5c8","Chloroflexi"="#84b59c","Cyanobacteria"="#427e83","Firmicutes"="#7997a1","Gemmatimonadetes"="#919f70","Nitrospirae"="#686a47",
                         "Planctomycetes"="#c8bd83","Proteobacteria"="#cb9e59","SPAM"="#ecc363","Tenericutes"="#c57b67","TM7"="#d65f3d","Verrucomicrobia"="#a04344","Other"="#bf7a8f")


rm$Phylum <- factor(rm$Phylum,levels=c("Actinobacteria","Chloroflexi","Cyanobacteria","Firmicutes","Bacteroidetes","Proteobacteria",
                                       "Acidobacteria","Armatimonadetes","Gemmatimonadetes","Nitrospirae",
                                       "Planctomycetes","SPAM","Tenericutes","TM7","Verrucomicrobia","Other" ))

data=subset(rm,Sample!="Post_flowering")
p<-ggplot(data, aes(x=Timepoint, y=Abundance, fill=Phylum)) + 
  geom_bar(stat = "identity", color = "NA")+facet_grid(SampleType~Sample,scales = "free_x",space="free_x")+ 
  scale_fill_manual(values = col_SampleByTreatment) +
  theme(axis.text.x=element_text(size=11,color="black",angle=90), 
        axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=11,face="bold"),
        text=element_text(size=11,face="bold"))+
  ylab("Relative Abundance")+xlab("")

ggsave(filename = "shotgun_whole_all-bacteria-genes-composition.jpg",plot = p,width=10,height=6)


##################figure2a-2b whole shotgun dots figure
setwd("~/Downloads/shotgun_whole_profile/hyper/")
library(ggplot2)
data<-read.table("whole_shotgun_rhi_soil_all_ratio_figure.txt",header=T,sep="\t")

ratio <- data[,c(5,7,8)]

ratio$Tissue<- factor(ratio$Tissue, levels=c("Rhi_up","Soil_up",
                                             "Rhi_down","Soil_down"))
#write.csv(rm.m,"pvalue.csv")
###
dat <- ratio[order(ratio[,3]),]
new_data = dat[order(dat$Ratio),]

new_data$Annotation <- factor(new_data$Annotation, 
                              levels=c("Mobilome_prophages_transposons",
                                       "Transcription",
                                       "Carbohydrate_transport_and_metabolism",
                                       "Secondary_metabolites_biosynthesis_transport_and_catabolism",
                                       "Inorganic_ion_transport_and_metabolism",
                                       "General_function_prediction_only",
                                       "Function_unknown",
                                       "Lipid_transport_and_metabolism",
                                       "Energy_production_and_conversion",
                                       "Replication_recombination_and_repair",
                                       "Amino_acid_transport_and_metabolism",
                                       "Nucleotide_transport_and_metabolism",
                                       "Coenzyme_transport_and_metabolism",
                                       "RNA_processing_and_modification",
                                       "Defense_mechanisms",
                                       "Cell_cycle_control_cell_division_chromosome_partitioning",
                                       "Signal_transduction_mechanisms",
                                       "Cell_wall/membrane/envelope_biogenesis",
                                       "Posttranslational_modification_protein_turnover_chaperones",
                                       "Translation_ribosomal_structure_and_biogenesis",
                                       "Cell_motility",
                                       "Intracellular_trafficking_secretion_and_vesicular_transport",
                                       "Extracellular_structures"))


p<-ggplot(new_data, aes(x=Ratio, y=Annotation,colour=Tissue)) +
  scale_colour_manual(name="",values = c("#004445","#6FB98F","#9B4F0F","#C99E10"))+
  scale_shape_manual(name="",values=c(1,1,19,19))+facet_wrap(~Tissue,ncol=2)+
  geom_point(size = 4,stroke = 1)+theme_bw()+
  geom_vline(aes(xintercept=1),colour="Black",size=1,linetype="dashed")+xlab("Enrichment ratio")+
  theme(axis.text.x=element_text(size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=10))+
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(strip.text = element_text(size=11,face="bold"))+ylab("")

ggsave(filename = "whole-shotgun-ratio_dots_figure.pdf",plot = p,width=15,height=8)
#####################################

