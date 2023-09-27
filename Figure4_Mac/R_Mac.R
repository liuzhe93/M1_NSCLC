setwd("D:/MyProjects/scRNA_immune/Mac_Lung/Figure4_Mac")
rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
packageVersion("Seurat")

load("D:/MyProjects/scRNA_immune/Mac_Lung/Figure3_scRNA/NSCLC.combined.pca14.res0.6.afteranno.RData")
levels(NSCLC.combined)
#[1] "Epithelial_cells"  "Monocyte"          "NK_cell"           "Macrophage"        "B_cell"            "Fibroblasts"      
#[7] "Endothelial_cells"      
Macrophage.subset<-subset(NSCLC.combined@meta.data, celltype=="Macrophage")
scRNAsub.Macrophage <- subset(NSCLC.combined, cells=row.names(Macrophage.subset))
scRNAsub.Macrophage
#An object of class Seurat 
#11994 features across 980 samples within 2 assays 
#Active assay: integrated (2000 features, 2000 variable features)
#1 other assay present: RNA
#3 dimensional reductions calculated: pca, umap, tsne
DefaultAssay(scRNAsub.Macrophage)<-'RNA'
scRNAsub.Macrophage
#An object of class Seurat 
#11994 features across 980 samples within 2 assays 
#Active assay: RNA (9994 features, 2000 variable features)
#1 other assay present: integrated
#3 dimensional reductions calculated: pca, umap, tsne
scRNAsub.Macrophage <- FindVariableFeatures(scRNAsub.Macrophage, selection.method = "vst", nfeatures = 2000)
scale.genes.Macrophage <-  rownames(scRNAsub.Macrophage)
scRNAsub.Macrophage <- ScaleData(scRNAsub.Macrophage)
scRNAsub.Macrophage <- RunPCA(scRNAsub.Macrophage, features = VariableFeatures(scRNAsub.Macrophage))

################################### determine the number of pcs##############################################
#筛选标准
#1.主成分累积贡献大于90%
#2.PC本身对方差贡献小于5%
#3.两个连续PCs之间差异小于0.1%
#将要检测的seurat对象传递给sce
library("ggplot2")
sce=scRNAsub.Macrophage
# Determine percent of variation associated with each PC
pct <- sce [["pca"]]@stdev / sum( sce [["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
pcs <- min(co1, co2)
pcs
#11
# Create a dataframe with values
plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
# Elbow plot to visualize 
pdf("Determine_NumberOfPCs.pdf")
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()
save(scRNAsub.Macrophage, file="Mac.combined_beforepcs.RData")

################################ determine the resolution###########################################################
library(clustree)
library(SeuratData)
scRNAsub.Macrophage <- FindVariableFeatures(scRNAsub.Macrophage, selection.method = "vst", nfeatures = 2000)
scRNAsub.Macrophage <- ScaleData(scRNAsub.Macrophage , verbose = FALSE)
scRNAsub.Macrophage <- RunPCA(scRNAsub.Macrophage)
scRNAsub.Macrophage <- RunUMAP(scRNAsub.Macrophage, reduction = "pca", dims = 1:11)
scRNAsub.Macrophage <- RunTSNE(scRNAsub.Macrophage, reduction = "pca", dims = 1:11)
scRNAsub.Macrophage <- FindNeighbors(scRNAsub.Macrophage, reduction = "pca", dims = 1:11)
scRNAsub.Macrophage <- FindClusters(
  object = scRNAsub.Macrophage,
  resolution = c(seq(0,1.6,.2))
)
clustree(scRNAsub.Macrophage@meta.data, prefix = "RNA_snn_res.")
pdf("Determine_resolution.pdf")
clustree(scRNAsub.Macrophage, prefix = "RNA_snn_res.") +
  guides(edge_colour = FALSE, edge_alpha = FALSE) +scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "blue", high = "red")+
  theme(legend.position = "bottom")
dev.off()
r=1.6
scRNAsub.Macrophage <- FindClusters(scRNAsub.Macrophage, resolution = r)
levels(scRNAsub.Macrophage)
#[1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11"
save(scRNAsub.Macrophage,file="Mac.res1.6.RData")
Mac.pca11.markers <- FindAllMarkers(object = scRNAsub.Macrophage, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                                    test.use = "wilcox")
save(Mac.pca11.markers,file="Mac.pca11.markers.res1.6.RData")
library(dplyr)
top30<-Mac.pca11.markers %>% group_by(cluster) %>% top_n(n=30,wt=avg_log2FC)
write.table(top30,"top30.pca11.res1.6.markers.csv",sep=",",quote=F)
head(Idents(scRNAsub.Macrophage), 5)
write.table(Mac.pca11.markers,"Mac.pca11.markers.csv",sep=",",quote=F)

# calculate the percentage and count the numbers
prop.table(table(Idents(scRNAsub.Macrophage), scRNAsub.Macrophage$label))
allsampleprop.each <-as.data.frame(prop.table(table(Idents(scRNAsub.Macrophage), scRNAsub.Macrophage$label)))
prop.table(table(Idents(scRNAsub.Macrophage)))
allsampleprop.total <-as.data.frame(prop.table(table(Idents(scRNAsub.Macrophage))))
write.csv(x = allsampleprop.each,file = 'anno.allsample.each.prop.csv',quote = T,row.names = T)
write.csv(x = allsampleprop.total,file = 'anno.allsample.total.prop.csv',quote = T,row.names = T)
table(Idents(scRNAsub.Macrophage))
#  0   1   2   3   4   5   6   7   8   9  10  11 
#124 123 117 107 106  96  66  57  54  49  45  36

require(org.Hs.eg.db)
library(topGO)
library(DOSE)
#devtools::install_github("eliocamp/ggnewscale")
library("ggnewscale")
x=as.list(org.Hs.egALIAS2EG)
geneList<-rep(0,nrow(scRNAsub.Macrophage))
names(geneList)<-row.names(scRNAsub.Macrophage)
geneList<-geneList[intersect(names(geneList),names(x))]
newwallgenes=names(geneList)
for (ii in 1:length(geneList)){
  names(geneList)[ii]<-x[[names(geneList)[ii]]][1]
  
}
gene_erichment_results=list()
for (c1 in as.character(unique(levels(Mac.pca11.markers$cluster)))){
  print(paste0("RUN ", c1))
  testgenes<-subset(Mac.pca11.markers,cluster==c1)$gene
  gene_erichment_results[[c1]]=list()
  testgeneList=geneList
  testgeneList[which(newwallgenes %in% testgenes)]= 1
  #gene_erichment_results=list()
  tab1=c()
  for(ont in c("BP","MF")){
    sampleGOdata<-suppressMessages(new("topGOdata",description="Simple session",ontology=ont,allGenes=as.factor(testgeneList),
                                       nodeSize=10,annot=annFUN.org,mapping="org.Hs.eg.db",ID="entrez"))
    resultTopGO.elim<-suppressMessages(runTest(sampleGOdata,algorithm="elim",statistic="Fisher"))
    
    resultTopGO.classic<-suppressMessages(runTest(sampleGOdata,algorithm="classic",statistic="Fisher"))
    tab1<-rbind(tab1,GenTable(sampleGOdata,Fisher.elim=resultTopGO.elim,Fisher.classic=resultTopGO.classic,orderBy="Fisher.elim",
                              topNodes=200))
  }
  gene_erichment_results[[c1]][["topGO"]]=tab1
  x<-suppressMessages(enrichDO(gene=names(testgeneList)[testgeneList==1],ont="DO",pvalueCutoff=1,pAdjustMethod="BH",universe=names(testgeneList),
                               minGSSize=5,maxGSSize=500,readable=T))
  gene_erichment_results[[c1]][["DO"]]=x
  dgn<-suppressMessages(enrichDGN(names(testgeneList)[testgeneList==1]))
  gene_erichment_results[[c1]][["DGN"]]=dgn
}
save(gene_erichment_results,file="Macrophage_gene_erichment_results.RData")
write.csv(gene_erichment_results[["0"]][["topGO"]],"Macrophage.cluster0.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["1"]][["topGO"]],"Macrophage.cluster1.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["2"]][["topGO"]],"Macrophage.cluster2.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["3"]][["topGO"]],"Macrophage.cluster3.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["4"]][["topGO"]],"Macrophage.cluster4.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["5"]][["topGO"]],"Macrophage.cluster5.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["6"]][["topGO"]],"Macrophage.cluster6.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["7"]][["topGO"]],"Macrophage.cluster7.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["8"]][["topGO"]],"Macrophage.cluster8.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["9"]][["topGO"]],"Macrophage.cluster9.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["10"]][["topGO"]],"Macrophage.cluster10.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["11"]][["topGO"]],"Macrophage.cluster11.GO.csv",quote=F,row.names=F)

#根据GO功能富集分析（促进免疫或者抑制免疫）+marker基因来注释macrophage亚群
new.cluster.ids<-c("M2", "M2", "M2", "M2", "M2", "M0", "M0", "M1", "M2", "M1","M2","M1")
names(new.cluster.ids) <- levels(scRNAsub.Macrophage)
scRNAsub.Macrophage <- RenameIdents(scRNAsub.Macrophage, new.cluster.ids)
scRNAsub.Macrophage$celltype<-Idents(scRNAsub.Macrophage)
scRNAsub.Macrophage
#An object of class Seurat 
#11994 features across 980 samples within 2 assays 
#Active assay: RNA (9994 features, 2000 variable features)
#1 other assay present: integrated
#3 dimensional reductions calculated: pca, umap, tsne

scRNAsub.Macrophage <- RunUMAP(scRNAsub.Macrophage, dims = 1:11)
scRNAsub.Macrophage <- RunTSNE(scRNAsub.Macrophage, dims = 1:11)
save(scRNAsub.Macrophage,file="scRNAsub.Macrophage.afteranno.RData")
scRNAsub.Macrophage.markers <- FindAllMarkers(object = scRNAsub.Macrophage, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(scRNAsub.Macrophage.markers,"scRNAsub.Macrophage.anno.csv",sep=",",quote=F)
save(scRNAsub.Macrophage.markers,file="scRNAsub.Macrophage.markers.afteranno.RData")
top30<-scRNAsub.Macrophage.markers %>% group_by(cluster) %>% top_n(n=30,wt=avg_log2FC)
write.table(top30,"scRNAsub.Macrophage.top30.anno.csv",sep=",",quote=F)
save(scRNAsub.Macrophage.markers,file="scRNAsub.Macrophage.markers.afteranno.RData")
# calculate the percentage and count the numbers
prop.table(table(Idents(scRNAsub.Macrophage), scRNAsub.Macrophage$label))
allsampleprop.each <-as.data.frame(prop.table(table(Idents(scRNAsub.Macrophage), scRNAsub.Macrophage$label)))
prop.table(table(Idents(scRNAsub.Macrophage)))
allsampleprop.total <-as.data.frame(prop.table(table(Idents(scRNAsub.Macrophage))))
write.csv(x = allsampleprop.each,file = 'annotation/anno.allsample.each.prop.csv',quote = T,row.names = T)
write.csv(x = allsampleprop.total,file = 'annotation/anno.allsample.total.prop.csv',quote = T,row.names = T)
table(Idents(scRNAsub.Macrophage))
# M2  M0  M1 
#676 162 142
pro.total <- table(Idents(scRNAsub.Macrophage),scRNAsub.Macrophage$label)
table(Idents(scRNAsub.Macrophage),scRNAsub.Macrophage$label)
pro.each <- table(Idents(scRNAsub.Macrophage),scRNAsub.Macrophage$label)
write.csv(x =pro.total,file = 'annotation/anno.pro.total.csv',quote = T,row.names = T)
write.csv(x =pro.each,file = 'annotation/anno.pro.each.csv',quote = T,row.names = T)
# visualization
scRNAsub.Macrophage$celltype <- Idents(scRNAsub.Macrophage)
pdf(paste0("annotation/tsne.pca11.res",r,".splitbyLabel_NoLegend.pdf"),width=40,height=10)
DimPlot(scRNAsub.Macrophage, reduction = "tsne", label = TRUE, pt.size=1, label.size = 0, split.by = 'label', group.by = "celltype")+
  NoLegend()
dev.off()
pdf(paste0("annotation/tsne.pca11.res",r,".splitbyLabel_Legend.pdf"),width=40,height=10)
DimPlot(scRNAsub.Macrophage, reduction = "tsne", label = TRUE, pt.size=1,label.size = 8, split.by = 'label', group.by = 'celltype')
dev.off()
pdf(paste0("annotation/tsne.pca11.res",r,"_Legend.pdf"),width=10,height=10)
DimPlot(scRNAsub.Macrophage, reduction = "tsne", label = TRUE, pt.size=1,label.size = 8, group.by = 'celltype')
dev.off()
pdf(paste0("annotation/tsne.pca11.res",r,"_NoLegend.pdf"),width=10,height=10)
DimPlot(scRNAsub.Macrophage, reduction = "tsne", label = FALSE, pt.size=1,label.size = 8, group.by = 'celltype')+NoLegend()
dev.off()
pdf("annotation/tsne.merged.pdf",width=10,height=10)
DimPlot(scRNAsub.Macrophage, reduction = "tsne", label = TRUE, pt.size=1, label.size = 8, group.by = "label")+NoLegend()
dev.off()

#########################################热图可视化##################################################################
pdf("heatmap.top30.pdf",width=24,height=18)
DoHeatmap(scRNAsub.Macrophage,features=top30$gene,cells = 1:500, size = 4, angle = 90, disp.min=-2, disp.max=2) + scale_fill_gradientn(colours=c("blue","white","red"))
dev.off()

#########################################气泡图可视化#################################################################
#从panglaodb数据库中挑选相应亚群的marker基因
features.plot <- c("MRC1", "PPARG", "TREM2", "FCGR3A", "FABP4", "CD52", 
                   "EIF3E", "TMEM176B", "TMEM176A", "VCAN", "SLC25A37", "PMP22",
                   "CD40", "CCL3", "CD86", "IL1B", "CSF1R", "IL6R")   

pdf("Mac.dittoDotPlot.pdf",width = 11, height = 6)
DotPlot(object = scRNAsub.Macrophage, features=features.plot,dot.scale = 6,cols = c("lightgrey", "red")) + RotatedAxis()
dev.off()

#####################################tsne可视化图各个亚群#############################################################
levels(scRNAsub.Macrophage)
#[1] "M2" "M0" "M1"
pdf("subpopulations/M2_tsne.pdf", width = 10, height = 10)
DimPlot(scRNAsub.Macrophage, reduction = "tsne",group.by ="celltype",cols=c('M2'='red','M0'='grey','M1'='grey'))+ NoLegend()
dev.off()

pdf("subpopulations/M0_tsne.pdf", width = 10, height = 10)
DimPlot(scRNAsub.Macrophage, reduction = "tsne",group.by ="celltype",cols=c('M2'='grey','M0'='red','M1'='grey'))+ NoLegend()
dev.off()

pdf("subpopulations/M1_tsne.pdf", width = 10, height = 10)
DimPlot(scRNAsub.Macrophage, reduction = "tsne",group.by ="celltype",cols=c('M2'='grey','M0'='grey','M1'='red'))+ NoLegend()
dev.off()

###############################堆积百分比柱形图#######################################################################
percentage<-read.csv("percentage/percentage.csv")
library("ggplot2")
library("reshape2")
mydata <- melt(percentage,id.vars="Celltype",variable.name="Label",value.name="CellNum")
pdf("percentage/percentage.pdf")
ggplot(mydata,aes(Celltype,CellNum,fill=Label))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("The fraction of cells")+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  coord_flip()
dev.off()
percentage_t<-read.csv("percentage/percentage_t.csv")
library("ggplot2")
library("reshape2")
mydata <- melt(percentage_t,id.vars="Celltype",variable.name="Label",value.name="CellNum")
pdf("percentage/percentage_t.pdf")
ggplot(mydata,aes(Celltype,CellNum,fill=Label))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("The fraction of cells")+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  coord_flip()
dev.off()

#############################单细胞功能注释和富集分析##############################################################
library("clusterProfiler")
library("org.Hs.eg.db")
ids=bitr(scRNAsub.Macrophage.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
sce.markers=merge(scRNAsub.Macrophage.markers,ids,by.x='gene',by.y='SYMBOL')
#View(sce.markers)
## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
## 它的返回值是一个列表，代表分组变量每个水平的观测。
gcSample=split(sce.markers$ENTREZID, sce.markers$cluster) 
#################################进行KEGG功能注释##################################################################
## KEGG
pdf("subpopulations/KEGG_enrichment.pdf", width = 6, height = 10)
xx <- compareCluster(gcSample,
                     fun = "enrichKEGG",
                     organism = "hsa", pvalueCutoff = 0.05
)
p <- dotplot(xx)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
dev.off()
ck <- setReadable(xx, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(ck) 
write.csv(ck,"subpopulations/enrichKEGG.csv")
###########################进行GO功能注释###############################################
## GO_BP
pdf("subpopulations/GOBP_enrichment.pdf", width = 6, height = 10)
xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05
)
p <- dotplot(xx)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
dev.off()
ck <- setReadable(xx, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(ck) 
write.csv(ck,"subpopulations/enrichGO_BP.csv")
## GO_MF
pdf("subpopulations/GOMF_enrichment.pdf", width = 6, height = 10)
xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05
)
p <- dotplot(xx)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
dev.off()
ck <- setReadable(xx, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(ck) 
write.csv(ck,"subpopulations/enrichGO_MF.csv")
## GO_CC
pdf("subpopulations/GOCC_enrichment.pdf", width = 6, height = 10)
xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05
)
p <- dotplot(xx)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
dev.off()
ck <- setReadable(xx, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(ck) 
write.csv(ck,"subpopulations/enrichGO_CC.csv")

#######################################优化图#################################################
setwd("D:/MyProjects/scRNA_immune/scRNA_CAF_PRAD/01_scRNA/subpopulations/EnrichmentAnalysis/")
df_GO<-read.csv("enrichGO_BP.csv", header = T)
#library("forcat")
library("tidyverse")

df_GO$Description<-as.factor(df_GO$Description)
df_GO$Description<-fct_inorder(df_GO$Description)

ggplot(df_GO, aes(Cluster, Description)) + 
  geom_point(aes(fill = p.adjust, size = Count), shape = 21) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text = element_text(color = "black", size = 10)) +
  scale_fill_gradient(low = "purple", high = "yellow") +
  labs(x = NULL, y = NULL) +
  coord_flip()

