#fig1b_jv
data=read.table("proteomics_1B_fold_log2.csv", row.names=1, header=TRUE, sep=",")

annotdf <- data.frame(category = c(rep("Complex II (SDH)", 4), rep("Complex III", 8), rep("Cytochrome C", 1), rep("Complex IV", 19), rep("Complex V", 15), rep("Tim22 Complex", 6), rep("Tim23 Complex", 4), rep("TOM Complex", 5), rep("ERMES", 5)))

rownames(annotdf) <- rownames(data)

library(pheatmap)
library(RColorBrewer)
pheatmap(as.matrix(t(data[,1:6])), scale='none', margins=c(75,1), col=rev(brewer.pal(11, "RdBu")), cluster_rows=FALSE, cluster_cols=FALSE, cellwidth=12, cellheight=40, gaps_row=c(3,6), gaps_col=c(4,12,13,32,47,53,57,62,67), breaks = c(-4.5,-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5,4.5), annotation_col = annotdf)

#fig2_jv
data2=read.table("proteomics_2.csv", header=TRUE, sep=",")
data2 = data2[1:47,]
rownames(data2) = data2[,1]
annotdf2 <- data.frame(category = c(rep("Complex II (SDH)", 4), rep("Complex III", 8), rep("Cytochrome C", 1), rep("Complex IV", 19), rep("Complex V", 15)))
rownames(annotdf2) <- rownames(data2)
pheatmap(as.matrix(t(data2[,2:10])), scale='none', margins=c(75,1), col=rev(brewer.pal(11, "RdBu")), cluster_rows=FALSE, cluster_cols=FALSE, cellwidth=12, cellheight=40, gaps_row=c(3,6,9), gaps_col=c(4,12,13,32,47), breaks = c(-4.5,-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5,4.5), annotation_col = annotdf2)

#proteomics_acetyl_28jul17
data = read.csv("Acetyl-CoA_heatmap.csv", header=TRUE, row.names=1, sep=",")

library(pheatmap)
library(RColorBrewer)
pheatmap(as.matrix(t(data[,1:9])), scale='none', margins=c(75,1), col=rev(brewer.pal(11, "RdBu")), cluster_rows=FALSE, cluster_cols=FALSE, cellwidth=12, cellheight=40, breaks = c(-4.5,-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5,4.5))



#heatmap_tn
setwd("~/Desktop/")
data_tn=read.table("Heat map mitoTn.csv", header=TRUE, sep=",")
rownames(data_tn) = data_tn[,1]
data_tn = data_tn[,-1]
library(pheatmap)
library(RColorBrewer)
pheatmap(as.matrix(t(data_tn),scale='none', margins=c(75,1), col=rev(brewer.pal(11, "RdBu")), cluster_rows=FALSE, cluster_cols=FALSE, cellwidth=12, cellheight=40, gaps_row=c(3,6), breaks = c(0,.2,.4,.6,.8,1,1.2,1.4,1.6,1.8)))


#heatmap_tn (the one that worked)
setwd("~/Desktop/")
library(pheatmap)
library(RColorBrewer)
data=read.table("~/Desktop/Heat map mitoTn.csv",header=TRUE,sep=',')
rownames(data) = data[,1]
pheatmap(as.matrix(t(data[,2:10])), scale='none', margins=c(75,1), col=rev(brewer.pal(11, "RdBu")), cluster_rows=FALSE, cluster_cols=FALSE, cellwidth=12, cellheight=40, gaps_row=c(3,6,9), breaks = c(0,.2,.4,.6,.8,1,1.2,1.4,1.6,1.8))
