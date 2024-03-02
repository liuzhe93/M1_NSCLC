setwd("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ExprData")
rm(list = ls())

################read expression profile data and merge#######################################
exp_GSE13213<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/external/GSE13213/exprSet_GSE13213.csv",header=T)
exp_GSE13213$gene_name<-exp_GSE13213$GENE_SYMBOL
exp_GSE13213<-exp_GSE13213[,-c(1,2,3)]
exp_GSE13213<-subset(exp_GSE13213, select=c("gene_name",colnames(exp_GSE13213)[-118]))
#exp_GSE13213: 117 samples 19412 genes

exp_GSE31210<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/external/GSE31210/exprSet_GSE31210.csv",header=T)
exp_GSE31210$gene_name<-exp_GSE31210$X
exp_GSE31210<-subset(exp_GSE31210, select=c("gene_name",colnames(exp_GSE31210)[-248]))
exp_GSE31210<-exp_GSE31210[,-2]
#exp_GSE31210: 246 samples 20857 genes

exp_LUAD<-read.table("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/TCGA/LUAD/TCGA_LUAD_anno_fpkm.txt",sep="\t",header=T)
exp_LUAD$gene_name<-row.names(exp_LUAD)
exp_LUAD<-subset(exp_LUAD, select=c("gene_name",colnames(exp_LUAD)[-597]))
#exp_LUAD: 596 samples 59427 genes

exp_LUSC<-read.table("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/TCGA/LUSC/TCGA_LUSC_anno_fpkm.txt",sep="\t",header=T)
exp_LUSC$gene_name<-row.names(exp_LUSC)
exp_LUSC<-subset(exp_LUSC, select=c("gene_name",colnames(exp_LUSC)[-552]))
#exp_LUSC: 551 samples 59427 genes

exp_GEO<-merge(exp_GSE13213, exp_GSE31210, by = "gene_name")
exp_TCGA<-merge(exp_LUAD, exp_LUSC, by = "gene_name")
exp_merged<-merge(exp_GEO, exp_TCGA, by = "gene_name")
#library("dplyr")
#xxxx<-exp_merged %>% filter(gene_name=="TTL")
expr_dat <- aggregate(x = exp_merged[,2:ncol(exp_merged)],
                      by = list(exp_merged$gene_name),
                      FUN = mean)
dim(expr_dat)
#[1] 15596  1511
rownames(expr_dat)<-expr_dat$Group.1
expr_dat<-expr_dat[,-1]
write.csv(expr_dat, "merged_expr_dat_with_batch_effect.csv", quote = F)

##########limma remove batch effect###################
library(limma)
library(Glimma)
library(edgeR)

pca_plot = function(dddd, gggg){
  library("FactoMineR")
  library("factoextra")
  df.pca <- PCA(t(dddd), graph = FALSE)
  fviz_pca_ind(df.pca,
               geom.ind = "point",
               col.ind = gggg,
               addEllipses = TRUE,
               legend.title = "Groups")
}
batch = c(c(rep("GSE13213",117),rep("GSE31210",246),rep("LUAD",596),rep("LUSC",551)))
df = limma::removeBatchEffect(expr_dat,
                              batch = batch)
pdf("RemoveBatchEffect_before_Batch.pdf")
g_batch=factor(batch)
pca_plot(expr_dat,g_batch)
dev.off()
pdf("RemoveBatchEffect_before_Type.pdf")
g_batch=factor(c(rep("NSCLC",117),rep("NSCLC",226),rep("Control",20),rep("NSCLC",537),rep("Control",59),rep("NSCLC",502),rep("Control",49)))
pca_plot(expr_dat,g_batch)
dev.off()
pdf("RemoveBatchEffect_after_Batch.pdf")
g_batch=factor(batch)
pca_plot(df,g_batch)
dev.off()
pdf("RemoveBatchEffect_after_Type.pdf")
g_batch=factor(c(rep("NSCLC",117),rep("NSCLC",226),rep("Control",20),rep("NSCLC",537),rep("Control",59),rep("NSCLC",502),rep("Control",49)))
pca_plot(df,g_batch)
dev.off()
write.csv(df, "merged_expr_dat_remove_batch_effect.csv", quote = F)


