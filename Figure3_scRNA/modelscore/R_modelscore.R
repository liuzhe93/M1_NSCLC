setwd("D:/MyProjects/scRNA_immune/Mac_Lung/Figure3_scRNA/modelscore")
rm(list=ls())

library("Seurat")
library("tidyverse")
library("msigdbr")
library("AUCell")

load("D:/MyProjects/scRNA_immune/Mac_Lung/Figure4_Mac/scRNAsub.Macrophage.afteranno.RData")
levels(scRNAsub.Macrophage)

DefaultAssay(scRNAsub.Macrophage) <- "RNA"
M1_features <- list(c("RUNX3", "ADAM19", "LIMD2", "PRKCB", "CYTIP", "FNBP1", "ICAM3", "WIPF1", "TAP1", "PSMB9", "LAP3"))

scRNAsub.Macrophage <- AddModuleScore(scRNAsub.Macrophage,
                                 features = M1_features,
                                 ctrl = 100,
                                 names = "M1_features")
head(scRNAsub.Macrophage@meta.data)

colnames(scRNAsub.Macrophage@meta.data)[38]<-"M1_Score"

cells_rankings <- AUCell_buildRankings(scRNAsub.Macrophage@assays$RNA@data, splitByBlocks = TRUE)
names(M1_features)<-"M1"
cells_AUC <- AUCell_calcAUC(M1_features, cells_rankings, aucMaxRank = nrow(cells_rankings)*0.1)
head(cells_AUC)
AUCell_auc <- as.numeric(getAUC(cells_AUC)["M1",])
scRNAsub.Macrophage$AUCell <- AUCell_auc
head(scRNAsub.Macrophage@meta.data)

pdf("AUCell_M1Score.pdf")
VlnPlot(scRNAsub.Macrophage, features = "AUCell", pt.size = 0, group.by = "celltype")
dev.off()


library("ggrepel")
tsne<-data.frame(scRNAsub.Macrophage@meta.data, scRNAsub.Macrophage@reductions$tsne@cell.embeddings)
head(tsne)
cell_type_med<-tsne %>%
  group_by(celltype) %>%
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )

pdf("TSNE_M1Score.pdf")
ggplot(tsne, aes(tSNE_1, tSNE_2)) +
  geom_point(aes(colour = AUCell)) +
  scale_color_gradient2(low = "green", mid = "white", high =  "red") +
  ggrepel::geom_label_repel(aes(label = celltype), fontface = "bold", data = cell_type_med, 
                            point.adding = unit(0.5, "lines")) + theme_bw()
dev.off()
  





