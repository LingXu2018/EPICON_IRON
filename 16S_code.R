setwd("~/Desktop/Code_input")
library("phyloseq")
library("ggplot2")
library("scales")
library("grid")
theme_set(theme_bw())
library("DESeq2")
library("ape")
library("vegan")
library("data.table")
#library("RColorBrewer")

Bushman=readRDS("Bushman_EPICON_16S.rds")
sample_data(Bushman)$SampleType<-factor(sample_data(Bushman)$SampleType, levels=c("Soil","Rhizosphere","Root"))
sample_data(Bushman)$Treatment<-factor(sample_data(Bushman)$Treatment, levels=c("Control","Pre_flowering","Post_flowering"))
sample_data(Bushman)$Timepoint<-factor(sample_data(Bushman)$Timepoint, levels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17"))
sample_data(Bushman)$TreatmentByTimepoint<-factor(sample_data(Bushman)$TreatmentByTimepoint, levels=c("ControlTP1",	"ControlTP2",	"ControlTP3",	"ControlTP4",	"ControlTP5",	"ControlTP6",	"ControlTP7",	"ControlTP8",	"ControlTP9",	"ControlTP10",	"ControlTP11",	"ControlTP12",	"ControlTP13",	"ControlTP14",	"ControlTP15",	"ControlTP16",	"ControlTP17",	"Pre_floweringTP3",	"Pre_floweringTP4",	"Pre_floweringTP5",	"Pre_floweringTP6",	"Pre_floweringTP7",	"Pre_floweringTP8",	"Pre_floweringTP9",	"Pre_floweringTP10",	"Pre_floweringTP11",	"Pre_floweringTP12",	"Pre_floweringTP13",	"Pre_floweringTP14",	"Pre_floweringTP15",	"Pre_floweringTP16",	"Pre_floweringTP17",	"Post_floweringTP10",	"Post_floweringTP11",	"Post_floweringTP12",	"Post_floweringTP13",	"Post_floweringTP14",	"Post_floweringTP15",	"Post_floweringTP16",	"Post_floweringTP17",	"Post_floweringTP8",	"Post_floweringTP9"	))

rar = rarefy_even_depth(Bushman,sample.size=13000)
rar = prune_taxa(taxa_sums(rar)>=1, rar)

#Fig.1a-1b, supp.2a=2b Relative abundance for different phyla
options(scipen=200)
b1_phylum <- tax_glom(Bushman, taxrank="Phylum")

c = subset_samples(b1_phylum,Treatment=="Control")
x <- c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17")
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
  control <- rbind(control,physeqdf[order(physeqdf$Phylum),])
}

control <- data.frame(control)
pre1 <- subset(control,Timepoint=="TP1"|Timepoint=="TP2")
pre1$Treatment[pre1$Treatment == "Control"] <- "Pre_flowering"

pre2 <- {}
pre = subset_samples(b1_phylum,Treatment=="Pre_flowering")
y <- c("TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17")
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
  pre2 <- rbind(pre2,physeqdf[order(physeqdf$Phylum),])
}

add2 <- subset(control,Treatment=="Control")
add3 <- subset(add2,Timepoint=="TP1"|Timepoint=="TP2"|Timepoint=="TP3"|Timepoint=="TP4"|Timepoint=="TP5"|Timepoint=="TP6"|Timepoint=="TP7")
add3$Treatment[add3$Treatment == "Control"] <- "Post_flowering"

p <- {}
y <- c("TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17")
post = subset_samples(b1_phylum,Treatment=="Post_flowering")
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
rm <- rbind(control,pre1,pre2,add3,p)

sum <- tapply(rm$Abundance,rm$Phylum, FUN=sum)
sum <- data.frame(sum)
sum <- sort(sum$sum,decreasing = T)
list <- c("Acidobacteria","Actinobacteria","Bacteroidetes","Chloroflexi","Cyanobacteria",
          "Firmicutes","Gemmatimonadetes","Nitrospirae",
          "Planctomycetes","Proteobacteria","Tenericutes","TM7","Verrucomicrobia")
a <- rm[rm$Phylum %in% list,]
b <- rm[!(rm$Phylum %in% list),]
b$Phylum <- "Other"
rm <- rbind(a,b)
rm$Treatment<-factor(rm$Treatment,levels=c("Control","Pre_flowering","Post_flowering"))
rm$Timepoint<-factor(rm$Timepoint, levels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17"))

col_SampleByTreatment<-c("Acidobacteria"="#4d5b6a","Actinobacteria"="#6a90b4","Armatimonadetes"="#599ec4","Bacteroidetes"="#a1c5c8","Chloroflexi"="#84b59c","Cyanobacteria"="#427e83","Firmicutes"="#7997a1","Gemmatimonadetes"="#919f70","Nitrospirae"="#686a47",
                         "Planctomycetes"="#c8bd83","Proteobacteria"="#cb9e59","SPAM"="#ecc363","Tenericutes"="#c57b67","TM7"="#d65f3d","Verrucomicrobia"="#a04344","Other"="#bf7a8f")

rm$Phylum <- factor(rm$Phylum,levels=c("Actinobacteria","Chloroflexi","Cyanobacteria","Firmicutes","Bacteroidetes","Proteobacteria",
                                       "Acidobacteria","Armatimonadetes","Gemmatimonadetes","Nitrospirae",
                                       "Planctomycetes","SPAM","Tenericutes","TM7","Verrucomicrobia","Other" ))

sub<-subset(rm,Treatment!="Post_flowering")
sub<-subset(sub,Timepoint=="TP3"|Timepoint=="TP4"|Timepoint=="TP8"|Timepoint=="TP9"|Timepoint=="TP10"|Timepoint=="TP11")

plot <- ggplot(sub, aes(x=Timepoint, y=Abundance, fill=Phylum)) + 
  geom_bar(stat = "identity", color = "NA")+facet_wrap(~Sample+Treatment,ncol=2)+ scale_fill_manual(values = col_SampleByTreatment) +
  theme(axis.text.x=element_text(size=11,color="black",angle=90), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=14),text=element_text(size=14))+
  ylab("Relative Abundance")+xlab("")

ggsave(filename = "16S_relative_all_root_Pre_newcolor.jpg",plot = plot,width=10,height=10)


