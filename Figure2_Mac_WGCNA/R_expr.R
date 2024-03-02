setwd("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ExprData/")
remove(list=ls())

counts_LUAD<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/TCGA/LUAD/count.csv", row.names = 1)
dim(counts_LUAD)
#60660   598

#######################LUAD data prepare######################
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

##########gene level filter: edgeR logCPM>1###################
dim(counts_LUAD_merged)
#[1] 60660   596
counts_LUAD<-as.data.frame(counts_LUAD_merged)
library(limma)
library(Glimma)
library(edgeR)
x<-as.matrix(counts_LUAD)
Treatment <- factor(c(rep("LUAD",537),rep("Control",59)))
y <- DGEList(counts = x, group = Treatment)
keep <- rowSums(cpm(y, log=TRUE) >= 1) > 1  #Pre-filtering ，过滤低表达基因
table(keep)
#keep
#FALSE  TRUE 
#32254 28406 
y <- y[keep,, keep.lib.sizes=FALSE]
dim(y)
#[1] 28406   596
dat = data.frame(y$counts)
write.csv(dat, "TCGA_LUAD_filter_logCPM_counts.csv", quote = F)

##################################From counts to FPKMs###############################################
#FPKM = read counts / (mapped reads (Millions) * exon length)
# exon length: kb
library("GenomicFeatures")
txdb <- makeTxDbFromGFF("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/TCGA/gencode.v36.annotation.gtf")
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene, function(x){sum(width(reduce(x)))})
class(exons_gene_lens)
length(exons_gene_lens)
exons_gene_lens <- as.data.frame(exons_gene_lens)
dim(exons_gene_lens)
exons_gene_lens1<-t(exons_gene_lens)
dim(exons_gene_lens1)
exons_gene_lens2 <- exons_gene_lens1[-1,]
exons_gene_lens2 <- as.data.frame(exons_gene_lens2)
rownames(exons_gene_lens2)
xc <- gsub("\\.(\\.?\\d*)", "", rownames(exons_gene_lens2))
rownames(exons_gene_lens2) = xc
colnames(exons_gene_lens2) = "Length"
exons_gene_lens2$X<-rownames(exons_gene_lens2)

counts <- read.csv("TCGA_LUAD_filter_logCPM_counts.csv", header = T)
counts$X <- gsub("\\.(\\.?\\d*)", "", counts$X)
count_with_length <- merge(counts, exons_gene_lens2, by = "X")
rownames(count_with_length) <- count_with_length$X
count_with_length <- count_with_length[,-1]
kb <- count_with_length$Length/1000
countdata <- count_with_length[, 1:596]
rpk <- countdata/kb
rpk
fpkm <- t(t(rpk)/colSums(countdata) * 10^6)
write.csv(fpkm, file = "TCGA_LUAD_count2fpkm.csv", quote = F)

############################expression profile gene annotation#############################################
library(rtracklayer)
x = rtracklayer::import("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/TCGA/gencode.v36.annotation.gtf.gz")
x2 = as.data.frame(x)
tj = as.data.frame(table(x2$type))
anno<-subset(x2,type == "gene", select = c("gene_name","gene_id"))
anno$gene_id <- gsub("\\.(\\.?\\d*)", "", anno$gene_id)
fpkm<-as.data.frame(fpkm)
fpkm$gene_id<-rownames(fpkm)
fpkm_anno<-merge(anno,fpkm,by="gene_id")
fpkm_anno<-fpkm_anno[,-1]
fpkm_anno<-aggregate(fpkm_anno,by=list(fpkm_anno$gene_name),FUN=median)
rownames(fpkm_anno)<-fpkm_anno$Group.1
fpkm_anno<-fpkm_anno[,c(-1,-2)]
dim(fpkm_anno)
#[1] 28308   596
write.table(fpkm_anno, file = "TCGA_LUAD_anno_count2fpkm.txt", sep = "\t", quote = F)

