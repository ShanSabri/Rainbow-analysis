##################################
## Setup
rm(list=ls(all=TRUE))
options(warn=-1)
setwd("~/your/working/directory")
source("utils.R")


##################################
## Read data
raw.data      <- ReadData('~/Dropbox/ErnstLab/Analyses/Ngoc_Nguyen/2017-07-17/data/reza-2016-sc-count.txt.gz')
filtered.data <- FilterGenes(raw.data, is.expr=1, n.cells=2)
md.cell       <- BuildCellMetaData(raw.data)
md.cell       <- subset(md.cell, Timepoint != 'P10')
md.cell       <- subset(md.cell, Genes >= 2000)
md.cell$Timepoint <- ifelse(md.cell$Timepoint %in% c('E12'), 'E12.5', ifelse(md.cell$Timepoint %in% c('E9'), 'E9.5', 'P1'))
md.cell$Timepoint <- factor(md.cell$Timepoint, levels=c('E9.5', 'E12.5', 'P1'))
normalized.data   <- NormalizeData(filtered.data[names(filtered.data) %in% row.names(md.cell)])
scaled.data       <- ScaleData(normalized.data)


##################################
# TSNE
npc <- 30
p   <- 27
i   <- 40
pca <- RunPCA(sclaed.data, n=npc)
projected.loadings <- t(type.of.data) %*% as.matrix(pca$v)
tsne    <- RunTSNE(projected.loadings, perp=p, id=i, iter=2000, theta=0, dims=2)
md.cell <-  cbind(md.cell[names(type.of.data),], data.frame(TSNE_1=tsne$Y[,1], TSNE_2=tsne$Y[,2]))


##################################
# Plot embedding 
require(ggplot2)
ggplot(md.cell, aes(TSNE_1,TSNE_2)) +
  geom_point(aes(size=Genes, colour=Timepoint)) + 
  theme_bw() +
  ggsave(file=paste0(getwd(), "/TSNE1-TSNE2.png"), height=8, width=8)


##################################
# K-means
set.seed(42)
k = 4
km.cluster <- as.factor(kmeans(data.frame(X=md.cell$TSNE_1, Y=md.cell$TSNE_2), k)$cluster)
md.cell <- cbind(md.cell, km.cluster)


##################################
# Plot embedding by K-means
ggplot(md.cell, aes(TSNE_1, TSNE_2)) + 
  geom_point(aes(size=Genes, colour=km.cluster, shape=Timepoint)) +
  theme_bw()+
  ggsave(paste0(getwd(), "/TSNE1-TSNE2-km.png"), height=5, width=6)


##################################
# Differental Expression
one.v.all <- DiffExp(md.cell, data=type.of.data, out=paste0(getwd(), "/diffexp/cluster.v.all/DEG", cores=7))
p.val.cutoff = 0.01 ; avg.diff.cutoff = 0.5
one.v.all <- one.v.all[lapply(one.v.all,length)>0]
sig.genes.list <- lapply(one.v.all, function(x) {subset(x, p.val<=p.val.cutoff & avg.diff>=avg.diff.cutoff)[order(subset(x, p.val<=p.val.cutoff & avg.diff>=avg.diff.cutoff)$avg.diff, decreasing=TRUE),]})
sig.genes.list <- sig.genes.list[sapply(sig.genes.list, function(x) dim(x)[1]) > 0]
sig.genes.list <- lapply(sig.genes.list, function(x) {cbind(x, gene=row.names(x), cluster=paste0("C", sapply(strsplit(colnames(x)[3], ".", fixed=TRUE), function(x) (x[3]))))})
sig.genes.list <- lapply(sig.genes.list, setNames, nm=c("p.val", "avg.diff", "pct.cluster", "pct.others", "mean.exp.cluster", "mean.exp.others", "log.mean.exp.cluster", "log.mean.exp.others", "gene", "cluster"))
sig.genes.one.v.all = do.call(rbind, sig.genes.list)


##################################
# Heatmap of DEGs
require(pheatmap)
require(RColorBrewer)
require(plyr)
require(dplyr)
top.n <- 100
sig.genes.one.v.all %>% group_by(cluster) %>% top_n(top.n, avg.diff) -> top.genes
top.genes$cluster <- factor(top.genes$cluster, levels=c('C1','C2','C3','C4'))
top.genes <- top.genes[with(top.genes, order(cluster, -pct.cluster)),]
sig.genes.one.v.all.dge <- type.of.data[as.vector((top.genes$gene)), row.names(km)]
heatmap.annotation <- md.cell[names(sig.genes.one.v.all.dge),c(3,9)]
heatmap.annotation <- heatmap.annotation[with(heatmap.annotation, order(KM_Cluster, Timepoint)),]
pheatmap(MinMax(sig.genes.one.v.all.dge[,row.names(heatmap.annotation)], min=-0.5, max=2),
         cluster_rows=FALSE, cluster_cols=FALSE,
         show_rownames=F, show_colnames=TRUE,
         color=colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100),         
         annotation_col=heatmap.annotation, 
         gaps_col=cumsum(table(md.cell$KM_Cluster)),
         gaps_row=cumsum(table(top.genes$cluster)),
         fontsize_row = 7,
         fontsize_col = 6, 
         filename="heatmap.png", height=10, width=15)