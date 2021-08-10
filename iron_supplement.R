library("ggplot2")
library(ggbeeswarm)

setwd("~/Downloads/Iron_MS/")

#figure7A
data <- read.csv("iron.txt", sep = "\t", 
                 skip = 0, header = TRUE,comment.char = "", check.names = FALSE)
data=data.frame(data)

data$Compound <- factor(data$Compound,levels=c("Blank","Low","High"))

p=ggplot(data, aes(x=Compound, y=Value, fill=Treatment,color=Treatment)) + geom_quasirandom(alpha=0.2,shape=1,size=0.8)+
  geom_boxplot(alpha=0.4)+ylab("Bacteria amount") +scale_color_manual(values = c("#427e83","#ecc363"), na.value = "white") + 
  scale_fill_manual(values=c("#427e83","#ecc363"),na.value="white")+
  theme_classic()+ facet_wrap(~Strain,nrow=2)
ggplot(ggsave(filename = "strep_pseudo_amount_qPCR_treatment.pdf",plot = p,width=4,height=4))


#Figure7B
#Phenotype
library("ggplot2")
library(ggbeeswarm)

setwd("~/Downloads/Iron_MS/")

data <- read.csv("iron_treatment.txt", sep = "\t", 
                 skip = 0, header = TRUE,comment.char = "", check.names = FALSE)
data=data.frame(data)
data$Compound <- factor(data$Compound,levels=c("No iron","low iron","high iron"))
p3=ggplot(data, aes(x=Compound, y=Rootweight, fill=Strain))+
  geom_boxplot(alpha=0.4)+ylab("root fresh weight (mg)") + facet_wrap(~Treatment,scales = "free_y",nrow=2)+ 
  theme_classic()+scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
ggplot(ggsave(filename = "root_Treatment.pdf",plot = p3,width=5,height=4))
