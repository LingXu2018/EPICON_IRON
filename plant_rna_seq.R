#Figure5
library("pheatmap")

setwd("~/Downloads/Iron_MS/")
rm=read.table("new_iron_genes_expression.txt",header=T,row.names = 1)
rm$Treatment <- factor(rm$Treatment,levels=c("Drought","Re-watering"))
rm$Timepoint<-factor(rm$Timepoint,levels=c("TP3","TP4","TP5","TP6","TP7",
                                           "TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15",
                                           "TP16","TP17"))
rm$Genotype<-factor(rm$Genotype,levels = c("RT430","BT642"))

annotation_col = rm[,c("Timepoint", "Treatment", "Genotype")]
annotation_col$Category=rownames(annotation_col)
data=data.frame(t(rm[,4:231]))

rownames(annotation_col)=colnames(data)
annotation_col=annotation_col[,-4]
out=pheatmap(data, scale = "row", cluster_cols = F, breaks = c(-10, seq(-2.5,2.5, length.out = 99), 10),
             show_colnames = TRUE, annotation_col = annotation_col,
             fontsize_row = 4,
             main = "")
res <- data[c(out$tree_row[["order"]]),out$tree_col[["order"]]]