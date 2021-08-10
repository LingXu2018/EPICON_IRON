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
library(ggpubr)

sub=readRDS("zmTOM_rhizosphere.rds")


#figure6A
###################################################
#plot pcoa by treatment and sampleType individually
##################################################
col_SampleByTreatment<-c("mediumpurple1","cadetblue3","lightsalmon2","lightskyblue","lightpink","lightblue1","darkseagreen","mistyrose2",
                         "blue","darkorchid1","darkolivegreen2","red","cornsilk",
                         "yellow","slategray1","plum")
sub.prop<-transform_sample_counts(sub, function(otu) otu/sum(otu))

zmTOM1_z_control=subset_samples(sub,Genotype=="zmTOM1"&SampleType=="Rhizosphere"&Treatment=="Control")
sub.prop2<-transform_sample_counts(zmTOM1_z_control, function(otu) otu/sum(otu))
ord.nmds.wunifrac2 <- ordinate(sub.prop2, method="CAP",distance="bray",formula=~SampleType+Genotype+Treatment) 
p1=phyloseq::plot_ordination(sub.prop2,ord.nmds.wunifrac2, color="Mutants") + 
  geom_point(size=3)+ggtitle("zmTOM1_z_control,P=0.0005")+theme(legend.position="none")+xlab("")+ylab("")+#facet_wrap(SampleType+Genotype~Treatment,scales="free",nrow=2)+
  theme(axis.text.x=element_text(size=11,color="black",angle=90), 
        axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=8,face="bold"),text=element_text(size=14,face="bold"))+
  scale_color_manual(values = col_SampleByTreatment)+stat_ellipse(type = "norm", linetype = 2)#+stat_ellipse(type = "t")
ggplot(ggsave(filename = "pcoa_maize_zmTOM1_z_control.pdf",plot = p2,width=4,height=4))

TP_dist <- phyloseq::distance(sub.prop2, method = "bray")
TP_data <- data.frame(sample_data(sub.prop2))
test <- adonis(TP_dist ~ Mutants , data = TP_data,permutations = 9999) #+ SampleType*Genotype +SampleType*Replicate +Genotype*Replicate
print(test)
#wunifrac#P=0.0004
#unifrac 0.0012
#bray 0.0005


zmTOM1_z_drought=subset_samples(sub,Genotype=="zmTOM1"&SampleType=="Rhizosphere"&Treatment=="Drought")
sub.prop4<-transform_sample_counts(zmTOM1_z_drought, function(otu) otu/sum(otu))
ord.nmds.wunifrac4 <- ordinate(sub.prop4, method="CAP",distance="bray",formula=~Treatment+Genotype+SampleType) 
p2=phyloseq::plot_ordination(sub.prop4,ord.nmds.wunifrac4, color="Mutants") + 
  geom_point(size=3)+ggtitle("zmTOM1_z_drought,P=0.5373")+theme(legend.position="none")+xlab("")+ylab("")+#facet_wrap(SampleType+Genotype~Treatment,scales="free",nrow=2)+
  theme(axis.text.x=element_text(size=11,color="black",angle=90), 
        axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=8,face="bold"),text=element_text(size=14,face="bold"))+
  scale_color_manual(values = col_SampleByTreatment)+stat_ellipse(type = "norm", linetype = 2)#+stat_ellipse(type = "t")
ggplot(ggsave(filename = "pcoa_maize_zmTOM1_z_drought.pdf",plot = p4,width=4,height=4))

TP_dist <- phyloseq::distance(sub.prop4, method = "bray")
TP_data <- data.frame(sample_data(sub.prop4))
test <- adonis(TP_dist ~ Mutants , data = TP_data,permutations = 9999) #+ SampleType*Genotype +SampleType*Replicate +Genotype*Replicate
print(test)
#wunifrac#P=0.7351
#unifrac 0.3655
#bray 0.5373
rich = estimate_richness(zmYS1_r_control)
pairwise.wilcox.test(rich$Shannon, sample_data(zmYS1_r_control)$Mutants)

ggarrange(p1,p2,labels = c("A", "B"),ncol = 2, nrow = 1)

#figure6B
########################relative abundance---Phylum#####################
options(scipen=200)
b1_phylum <- tax_glom(sub, taxrank="Phylum")

c = subset_samples(b1_phylum,Treatment=="Control" & Genotype=="zmTOM1")
x <- c("Mutant","Wild")
control_zmTOM1 <- {}
for(i in x){
  TP = subset_samples(c, Mutants == i)
  merged_r = merge_samples(TP, "SampleType")
  merged_r_m = transform_sample_counts(merged_r , function(x) 100 * x/sum(x))
  physeqdf <- psmelt(merged_r_m)
  temp <- c("Rhizosphere","Root")
  physeqdf$Sample <- factor(physeqdf$Sample, levels = c(temp))
  n = nrow(physeqdf[order(physeqdf$Phylum),])
  x = sample(i,n,replace=TRUE)
  ST <- sample("Control",n,replace=TRUE)
  physeqdf[order(physeqdf$Phylum),]$Mutants <- x 
  physeqdf[order(physeqdf$Phylum),]$Treatment <- ST
  physeqdf[order(physeqdf$Phylum),]$Genotype <- "zmTOM1"
  physeqdf$Phylum <- as.character(physeqdf$Phylum)
  #physeqdf$Phylum[30:129] <- "Other"
  control_zmTOM1 <- rbind(control_zmTOM1,physeqdf[order(physeqdf$Phylum),])
}
control_zmTOM1 <- data.frame(control_zmTOM1)


d = subset_samples(b1_phylum,Treatment=="Drought" & Genotype=="zmTOM1")
x <- c("Mutant","Wild")
drought_zmTOM1 <- {}
for(i in x){
  TP = subset_samples(d, Mutants == i)
  merged_r = merge_samples(TP, "SampleType")
  merged_r_m = transform_sample_counts(merged_r , function(x) 100 * x/sum(x))
  physeqdf <- psmelt(merged_r_m)
  temp <- c("Rhizosphere","Root")
  physeqdf$Sample <- factor(physeqdf$Sample, levels = c(temp))
  n = nrow(physeqdf[order(physeqdf$Phylum),])
  x = sample(i,n,replace=TRUE)
  ST <- sample("Drought",n,replace=TRUE)
  physeqdf[order(physeqdf$Phylum),]$Mutants <- x 
  physeqdf[order(physeqdf$Phylum),]$Treatment <- ST
  physeqdf[order(physeqdf$Phylum),]$Genotype <- "zmTOM1"
  physeqdf$Phylum <- as.character(physeqdf$Phylum)
  drought_zmTOM1 <- rbind(drought_zmTOM1,physeqdf[order(physeqdf$Phylum),])
}
drought_zmTOM1 <- data.frame(drought_zmTOM1)


rm <- rbind(control_zmTOM1,drought_zmTOM1)

sum <- tapply(rm$Abundance,rm$Phylum, FUN=sum)
sum <- data.frame(sum)
sum <- sort(sum$sum,decreasing = T)
#list <- rownames(sum[1:10])
list <- c("Acidobacteria","Actinobacteria","Bacteroidetes","Chloroflexi","Cyanobacteria",
          "Firmicutes","Gemmatimonadetes","Nitrospirae",
          "Planctomycetes","Proteobacteria","Tenericutes","TM7","Verrucomicrobia")
a <- rm[rm$Phylum %in% list,]
b <- rm[!(rm$Phylum %in% list),]
b$Phylum <- "Other"
rm <- rbind(a,b)
rm$Treatment<-factor(rm$Treatment,levels=c("Control","Drought"))

rm$Phylum <- factor(rm$Phylum,levels=c("Actinobacteria","Chloroflexi","Cyanobacteria","Firmicutes","Bacteroidetes","Proteobacteria",
                                       "Acidobacteria","Armatimonadetes","Gemmatimonadetes","Nitrospirae",
                                       "Planctomycetes","SPAM","Tenericutes","TM7","Verrucomicrobia","Other" ))


plot <- ggplot(rm, aes(x=Mutants, y=Abundance, fill=Phylum)) + 
  geom_bar(stat = "identity", color = "NA")+facet_grid(Treatment~Sample+Genotype)+ scale_fill_manual(values = col_SampleByTreatment) +
  theme(axis.text.x=element_text(size=11,color="black",angle=90,face="bold"), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=14,face="bold"),text=element_text(size=14,face="bold"))+
  ylab("Relative Abundance")+xlab("")

ggsave(filename = "relative_all_maize_phylum.jpg",plot = plot,width=8,height=6)

#figure6C
actino=subset(rm,Phylum=="Actinobacteria")
actino$Mutants<-factor(actino$Mutants,levels=c("Wild","Mutant"))
actino_new=actino[,c("Sample","Abundance","Genotype","Mutants","Treatment")]

library(dplyr)
library(tidyr)

dat=actino_new %>%
  spread(Mutants, Abundance) %>%
  mutate(ratio = Mutant/Wild)

ggplot(dat, aes(x="", y=ratio, fill="#6a90b4")) + 
  geom_bar(stat = "identity", color = "NA")+facet_grid(Genotype~Sample+Treatment)+ 
  theme(axis.text.x=element_text(size=11,color="black",angle=90), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=14),text=element_text(size=14))+
  ylab("Ratio of Relative Abundance for Actinobacteria (mutant/wild)")+xlab("")+scale_fill_manual(values="#6a90b4")
