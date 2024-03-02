setwd("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ExprData/")
remove(list=ls())

fpkms_LUAD<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/TCGA/LUAD/FPKM.csv", row.names = 1)
dim(fpkms_LUAD)
#60660   598

#######################LUAD data prepare######################
fpkms_LUAD_t<-t(fpkms_LUAD)
fpkms_LUAD_t<-as.data.frame(fpkms_LUAD_t)
cancer_bar<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/TCGA/LUAD/Cancer_barcode.csv")
cancer_bar$x<-gsub("-","\\.",cancer_bar$x)
fpkms_LUAD_t_can<-fpkms_LUAD_t[cancer_bar$x,]
normal_bar<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/TCGA/LUAD/Normal_barcode.csv")
normal_bar$x<-gsub("-","\\.",normal_bar$x)
fpkms_LUAD_t_nor<-fpkms_LUAD_t[normal_bar$x,]
dim(fpkms_LUAD_t_can)
#[1]   537 60660
dim(fpkms_LUAD_t_nor)
#[1]    59 60660
# in total, there are 537 cancer samples and 59 normal samples.
fpkms_LUAD_t_merged<-rbind(fpkms_LUAD_t_can, fpkms_LUAD_t_nor)
fpkms_LUAD_merged<-t(fpkms_LUAD_t_merged)
dim(fpkms_LUAD_merged)
#[1] 60660   596
fpkms_LUAD<-as.data.frame(fpkms_LUAD_merged)
write.csv(fpkms_LUAD, file = "TCGA_LUAD_FPKM.csv", quote = F)

############################expression profile gene annotation#############################################
library(rtracklayer)
x = rtracklayer::import("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/TCGA/gencode.v36.annotation.gtf.gz")
x2 = as.data.frame(x)
tj = as.data.frame(table(x2$type))
anno<-subset(x2,type == "gene", select = c("gene_name","gene_id"))
fpkms_LUAD<-read.csv("TCGA_LUAD_FPKM.csv",header=T,row.names = 1)
fpkm<-as.data.frame(fpkms_LUAD)
fpkm$gene_id<-rownames(fpkm)
fpkm_anno<-merge(anno,fpkm,by="gene_id")
fpkm_anno<-fpkm_anno[,-1]
fpkm_anno<-aggregate(fpkm_anno,by=list(fpkm_anno$gene_name),FUN=median)
rownames(fpkm_anno)<-fpkm_anno$Group.1
fpkm_anno<-fpkm_anno[,c(-1,-2)]
dim(fpkm_anno)
#[1] 59427   596
write.table(fpkm_anno, file = "TCGA_LUAD_anno_fpkm.txt", sep = "\t", quote = F)

