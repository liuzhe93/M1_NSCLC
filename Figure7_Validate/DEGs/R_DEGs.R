setwd("D:/MyProjects/scRNA_immune/Mac_Lung/Figure7_Validate/DEGs")
rm(list = ls())

risk<-read.table("D:/MyProjects/scRNA_immune/Mac_Lung/Figure5_M1_prognosis/risk.txt",
                 sep="\t",header=T, row.names = 1)
head(risk)
risk$sample_name<-row.names(risk)

for(i in 1:nrow(risk)){
  tmp <- strsplit(risk$sample_name[i], split = "_")[[1]][1]
  risk$sample_name[i] <- tmp
}

risk<-subset(risk, select = c(risk, sample_name))

#######################LUAD data prepare######################
counts_LUAD<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/TCGA/LUAD/count.csv", row.names = 1)
dim(counts_LUAD)
#60660   598
counts_LUAD_t<-t(counts_LUAD)
counts_LUAD_t<-as.data.frame(counts_LUAD_t)
cancer_bar<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/TCGA/LUAD/Cancer_barcode.csv")
cancer_bar$x<-gsub("-","\\.",cancer_bar$x)
counts_LUAD_t_can<-counts_LUAD_t[cancer_bar$x,]
normal_bar<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/TCGA/LUAD/Normal_barcode.csv")
normal_bar$x<-gsub("-","\\.",normal_bar$x)
counts_LUAD_t_nor<-counts_LUAD_t[normal_bar$x,]
dim(counts_LUAD_t_can)
#[1]   537 60660
dim(counts_LUAD_t_nor)
#[1]    59 60660
# in total, there are 537 cancer samples and 59 normal samples.
counts_LUAD_t_merged<-rbind(counts_LUAD_t_can, counts_LUAD_t_nor)
counts_LUAD_merged<-t(counts_LUAD_t_merged)
counts_LUAD_merged<-as.data.frame(counts_LUAD_merged)
counts_LUAD_merged$X<-rownames(counts_LUAD_merged)
counts_LUAD_merged$X <- gsub("\\.(\\.?\\d*)", "", counts_LUAD_merged$X)

library(rtracklayer)
x = rtracklayer::import("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/TCGA/gencode.v36.annotation.gtf.gz")
x2 = as.data.frame(x)
tj = as.data.frame(table(x2$type))
anno<-subset(x2,type == "gene", select = c("gene_name","gene_id"))
anno$gene_id <- gsub("\\.(\\.?\\d*)", "", anno$gene_id)
LUAD_anno<-merge(anno, counts_LUAD_merged, by.x = "gene_id", by.y="X")
LUAD_anno<-LUAD_anno[,-1]
LUAD_anno<-aggregate(LUAD_anno,by=list(LUAD_anno$gene_name),FUN=median)
rownames(LUAD_anno)<-LUAD_anno$Group.1
LUAD_anno<-LUAD_anno[,c(-1,-2)]
dim(LUAD_anno)
#[1] 59427   596
write.table(LUAD_anno, file = "TCGA_LUAD_anno_count.txt", sep = "\t", quote = F)

#######################LUSC data prepare######################
counts_LUSC<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/TCGA/LUSC/count.csv", row.names = 1)
dim(counts_LUSC)
#60660   551
counts_LUSC_t<-t(counts_LUSC)
counts_LUSC_t<-as.data.frame(counts_LUSC_t)
cancer_bar<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/TCGA/LUSC/Cancer_barcode.csv")
cancer_bar$x<-gsub("-","\\.",cancer_bar$x)
counts_LUSC_t_can<-counts_LUSC_t[cancer_bar$x,]
normal_bar<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/TCGA/LUSC/Normal_barcode.csv")
normal_bar$x<-gsub("-","\\.",normal_bar$x)
counts_LUSC_t_nor<-counts_LUSC_t[normal_bar$x,]
dim(counts_LUSC_t_can)
#[1]   502 60660
dim(counts_LUSC_t_nor)
#[1]    49 60660
# in total, there are 502 cancer samples and 49 normal samples.
counts_LUSC_t_merged<-rbind(counts_LUSC_t_can, counts_LUSC_t_nor)
counts_LUSC_merged<-t(counts_LUSC_t_merged)
counts_LUSC_merged<-as.data.frame(counts_LUSC_merged)
counts_LUSC_merged$X<-rownames(counts_LUSC_merged)
counts_LUSC_merged$X <- gsub("\\.(\\.?\\d*)", "", counts_LUSC_merged$X)
LUSC_anno<-merge(anno, counts_LUSC_merged, by.x = "gene_id", by.y="X")
LUSC_anno<-LUSC_anno[,-1]
LUSC_anno<-aggregate(LUSC_anno,by=list(LUSC_anno$gene_name),FUN=median)
rownames(LUSC_anno)<-LUSC_anno$Group.1
LUSC_anno<-LUSC_anno[,c(-1,-2)]
dim(LUSC_anno)
#[1] 59427   551
write.table(LUSC_anno, file = "TCGA_LUSC_anno_count.txt", sep = "\t", quote = F)
#########################################

exp_TCGA<-cbind(LUAD_anno, LUSC_anno)
dataexpr_t<-t(exp_TCGA)
dataexpr_t<-as.data.frame(dataexpr_t)
dataexpr_t$sample_name<-row.names(dataexpr_t)
dataexpr_t$sample_name<-substr(dataexpr_t$sample_name,1,12)
dataexpr_t$sample_name<-gsub("\\.","-",dataexpr_t$sample_name)
merged_data<-merge(dataexpr_t, risk, by = "sample_name")
dim(dataexpr_t)
#[1]  1147 59428
dim(risk)
#[1] 1475    2
dim(merged_data)
#[1]  1430 59429

merged_data<-merged_data[order(merged_data$risk),]
merged_data[1:6,59426:59429]
table(merged_data$risk)
#high  low 
#665  765 

for(i in 1:1430){
  merged_data$sample_name[i]<-paste(merged_data$sample_name[i],i,sep = "_")
}
row.names(merged_data)<-merged_data$sample_name
merged_data<-merged_data[,-1]
merged_data<-merged_data[,-59428]
dim(merged_data)
#[1]  1430 59427
# the first 665 rows indicate the high-risk groups 
# the following 765 rows indicate the low-risk groups
merged_data<-t(merged_data)
write.csv(merged_data, "Count_highlow_risk.csv", quote = F)

Group <- factor(rep(c("high","low"),times=c(665,765)),levels = c("low","high"))
table(Group)
#Group
#low high 
#765  665
design <- model.matrix(~0+Group)
colnames(design)= levels(Group)
rownames(design)=colnames(merged_data)
write.csv(merged_data,"dataexpr.csv",quote = F)
#merged_data<-read.csv("dataexpr.csv",row.names=1)

library("limma")
library("tidyverse")
library("stringr")
library("edgeR")
table(is.na(merged_data))
#FALSE 
#84980610 
dge <- DGEList(counts=merged_data)
dge <- calcNormFactors(dge)

v <- voom(dge,design, normalize="quantile")
fit <- lmFit(v, design)
constrasts = paste(rev(levels(Group)),collapse = "-")
cont.matrix <- makeContrasts(contrasts=constrasts,levels = design) 
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)
DEG = topTable(fit2, coef=constrasts, n=Inf)
DEG = na.omit(DEG)
logFC_t=0.5849625
P.Value_t = 0.05
k1 = (DEG$P.Value < P.Value_t)&(DEG$logFC < -logFC_t)
k2 = (DEG$P.Value < P.Value_t)&(DEG$logFC > logFC_t)
change = ifelse(k1,"DOWN",ifelse(k2,"UP","stable"))
DEG$change <- change
DEG_limma <- DEG
table(DEG_limma$change)
#DOWN stable     UP 
# 239  59102     86
write.csv(DEG_limma,"deg_High_LowRisk.csv",quote=F)

up_genes<-subset(DEG_limma, change == "UP")
up_list<-rownames(up_genes)
write.table(up_list, "up_genes.txt", sep = "\t", row.names = F, col.names = F, quote = F)

down_genes<-subset(DEG_limma, change == "DOWN")
down_list<-rownames(down_genes)
write.table(down_list, "down_genes.txt", sep = "\t", row.names = F, col.names = F, quote = F)

library("ggplot2")
deg<-DEG_limma
head(deg)
deg$color<-ifelse(deg$P.Value<0.05 & abs(deg$logFC)>logFC_t, ifelse(deg$logFC< -logFC_t, "blue", "red"),"gray")
color<-c(red = "red", gray = "gray", blue = "blue")
p<-ggplot(deg, aes(logFC, -log10(P.Value), col = color)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = color) +
  labs(x="log2 (fold change)", y = "-log10(P-Value)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
pdf("DEG_High_Low_VP.pdf")
p
dev.off()


#https://metascape.org/gp/index.html#/main/step1
#Figure 1. Bar graph of enriched terms across input gene lists, colored by p-values.

