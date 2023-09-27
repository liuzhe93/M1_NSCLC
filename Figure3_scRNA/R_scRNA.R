setwd("/Users/liuzhe/Desktop/scRNA_immune/Mac_Lung/Figure3_scRNA")
rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
packageVersion("Seurat")
library("DoubletFinder")


##################第1个病人：GSM3304007_P1_Tumorz########################################
ScRNA_exp_P1 <- read.table(
  "/Users/liuzhe/Desktop/scRNA_immune/Mac_Lung//DataCollection/scRNA-seq/GSE117570_RAW/GSM3304007_P1_Tumor_processed_data.txt.gz",
  row.names = 1,
  header = T)
Seurat_object_P1 <- CreateSeuratObject(
  counts = ScRNA_exp_P1, 
  min.cells = 3, 
  min.features = 200)
Seurat_object_P1 <- NormalizeData(Seurat_object_P1)
Seurat_object_P1 <- FindVariableFeatures(Seurat_object_P1, nfeatures = 2000)
Seurat_object_P1[["percent.mt"]] <- PercentageFeatureSet(Seurat_object_P1, pattern = "^MT-")
Seurat_object_P1
Seurat_object_P1 <- ScaleData(Seurat_object_P1)
Seurat_object_P1 <- RunPCA(Seurat_object_P1)
Seurat_object_P1 <- RunUMAP(Seurat_object_P1, reduction = "pca", dims = 1:20)
Seurat_object_P1 <- RunTSNE(Seurat_object_P1, reduction = "pca", dims = 1:20)
## Pre-process Seurat object (sctransform) -----------------------------------------------------------------------------------
Seurat_object_P1 <- SCTransform(Seurat_object_P1)
Seurat_object_P1 <- RunPCA(Seurat_object_P1)
Seurat_object_P1 <- RunUMAP(Seurat_object_P1, reduction = "pca", dims=1:20)
Seurat_object_P1 <- RunTSNE(Seurat_object_P1, reduction = "pca", dims=1:20)
sweep.res.list_brain <- paramSweep_v3(Seurat_object_P1, PCs = 1:20, sct = FALSE)
sweep.stats_brain <- summarizeSweep(sweep.res.list_brain, GT = FALSE)
head(sweep.stats_brain)
sweep.stats_brain[order(sweep.stats_brain$BCreal),]
bcmvn_brain <- find.pK(sweep.stats_brain)
mpK<-as.numeric(as.vector(bcmvn_brain$pK[which.max(bcmvn_brain$BCmetric)]))
DoubletRate = ncol(Seurat_object_P1)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
#DoubletRate = 0.075    # 直接查表，10000细胞对应的doublets rate是～7.6%
#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
Seurat_object_P1 <- FindNeighbors(Seurat_object_P1, reduction = "pca", dims = 1:20)
Seurat_object_P1 <- FindClusters(Seurat_object_P1, resolution = 0.4)
levels(Seurat_object_P1)
#[1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11"
Seurat_object_P1 <- RunUMAP(Seurat_object_P1, reduction = "pca", dims = 1:20)
Seurat_object_P1 <- RunTSNE(Seurat_object_P1, reduction = "pca", dims = 1:20)
annotations <- Seurat_object_P1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(DoubletRate*nrow(Seurat_object_P1@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
Seurat_object_P1 <- doubletFinder_v3(Seurat_object_P1, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
Seurat_object_P1 <- doubletFinder_v3(Seurat_object_P1, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.14_27", sct = FALSE)
pdf("DoubletFinder_P1.pdf")
DimPlot(Seurat_object_P1, reduction = "tsne", group.by = "DF.classifications_0.25_0.14_27", pt.size = 0.6)
dev.off()
Seurat_object_P1@meta.data$singledouble<-Seurat_object_P1@meta.data$'DF.classifications_0.25_0.14_27'
Seurat_object_P1.singlet <- subset(Seurat_object_P1, subset = singledouble == "Singlet")
Seurat_object_P1$'DF.classifications_0.25_0.14_27'<-Idents(Seurat_object_P1)
DimPlot(Seurat_object_P1, reduction = "tsne", label = TRUE, pt.size=0.6,label.size = 8, group.by = "DF.classifications_0.25_0.14_27")
Seurat_object_P1_singlet<-Seurat_object_P1.singlet
plot(x=Seurat_object_P1_singlet@meta.data$nCount_RNA,y=Seurat_object_P1_singlet@meta.data$nFeature_RNA)
# check the metadata in the new Seurat objects
head(Seurat_object_P1_singlet@meta.data)
tail(Seurat_object_P1_singlet@meta.data)
# Create .RData object to load at any time
save(Seurat_object_P1_singlet, file="Seurat_object_P1_singlet.RData")
Seurat_object_P1_singlet$log10GenesPerUMI <- log10(Seurat_object_P1_singlet$nFeature_RNA) / log10(Seurat_object_P1_singlet$nCount_RNA)
Seurat_object_P1_singlet$mitoRatio <- PercentageFeatureSet(object = Seurat_object_P1_singlet, pattern = "^MT-")
Seurat_object_P1_singlet$mitoRatio <- Seurat_object_P1_singlet@meta.data$mitoRatio / 100
Seurat_object_P1_singletmetadata <- Seurat_object_P1_singlet@meta.data
Seurat_object_P1_singletmetadata$cells <- rownames(Seurat_object_P1_singletmetadata)
Seurat_object_P1_singletmetadata <- Seurat_object_P1_singletmetadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
Seurat_object_P1_singlet
Seurat_object_P1_singlet@meta.data <- Seurat_object_P1_singletmetadata
counts <- GetAssayData(object = Seurat_object_P1_singlet, slot = "counts")
Seurat_object_P1_singlet <- CreateSeuratObject(counts, meta.data = Seurat_object_P1_singlet@meta.data)
Seurat_object_P1_singlet$label <- "P1"
P1_norm <- NormalizeData(Seurat_object_P1_singlet, normalization.method = "LogNormalize", scale.factor = 10000)
P1_norm <- FindVariableFeatures(P1_norm, selection.method = "vst", nfeatures = 2000)
Seurat_object_P1_singlet$log10GenesPerUMI <- log10(Seurat_object_P1_singlet$nFeature_RNA) / log10(Seurat_object_P1_singlet$nCount_RNA)
Seurat_object_P1_singlet$mitoRatio <- PercentageFeatureSet(object = Seurat_object_P1_singlet, pattern = "^MT-")
Seurat_object_P1_singlet$mitoRatio <- Seurat_object_P1_singlet@meta.data$mitoRatio / 100
VlnPlot(Seurat_object_P1_singlet, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3)
plot1 <- FeatureScatter(Seurat_object_P1_singlet, feature1 = "nCount_RNA", feature2 = "mitoRatio")
plot2 <- FeatureScatter(Seurat_object_P1_singlet, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
filtered_P1 <- subset(x = Seurat_object_P1_singlet, 
                          subset= (nCount_RNA < 6000) & 
                            (nFeature_RNA > 100) &
                            (nFeature_RNA < 2500) & 
                            (log10GenesPerUMI > 0.70) & 
                            (mitoRatio < 0.3))
Seurat_object_P1
#An object of class Seurat 
#6611 features across 1832 samples within 2 assays 
#Active assay: SCT (6611 features, 3000 variable features)
#1 other assay present: RNA
#3 dimensional reductions calculated: pca, umap, tsne
Seurat_object_P1_singlet
#An object of class Seurat 
#6611 features across 1805 samples within 1 assay 
#Active assay: RNA (6611 features, 0 variable features)
filtered_P1
#An object of class Seurat 
#6611 features across 1788 samples within 1 assay 
#Active assay: RNA (6611 features, 0 variable features)
counts <- GetAssayData(object = filtered_P1, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 20
filtered_counts <- counts[keep_genes, ]
filtered_P1 <- CreateSeuratObject(filtered_counts, meta.data = Seurat_object_P1_singlet@meta.data)
filtered_P1$label <- "P1"
filtered_P1
#An object of class Seurat 
#6610 features across 1788 samples within 1 assay 
#Active assay: RNA (6610 features, 0 variable features)
save(filtered_P1, file="filtered_P1.RData")
filtered_P1_norm<-NormalizeData(filtered_P1)

##################第2个病人：GSM3304009_P2_Tumor#########################################
ScRNA_exp_P2 <- read.table(
  "D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/scRNA-seq/GSE117570_RAW/GSM3304009_P2_Tumor_processed_data.txt.gz",
  row.names = 1,
  header = T)
Seurat_object_P2 <- CreateSeuratObject(
  counts = ScRNA_exp_P2, 
  min.cells = 3, 
  min.features = 200)
Seurat_object_P2 <- NormalizeData(Seurat_object_P2)
Seurat_object_P2 <- FindVariableFeatures(Seurat_object_P2, nfeatures = 2000)
Seurat_object_P2[["percent.mt"]] <- PercentageFeatureSet(Seurat_object_P2, pattern = "^MT-")
Seurat_object_P2
Seurat_object_P2 <- ScaleData(Seurat_object_P2)
Seurat_object_P2 <- RunPCA(Seurat_object_P2)
Seurat_object_P2 <- RunUMAP(Seurat_object_P2, reduction = "pca", dims = 1:20)
Seurat_object_P2 <- RunTSNE(Seurat_object_P2, reduction = "pca", dims = 1:20)
## Pre-process Seurat object (sctransform) -----------------------------------------------------------------------------------
Seurat_object_P2 <- SCTransform(Seurat_object_P2)
Seurat_object_P2 <- RunPCA(Seurat_object_P2)
Seurat_object_P2 <- RunUMAP(Seurat_object_P2, reduction = "pca", dims=1:20)
Seurat_object_P2 <- RunTSNE(Seurat_object_P2, reduction = "pca", dims=1:20)
sweep.res.list_brain <- paramSweep_v3(Seurat_object_P2, PCs = 1:20, sct = FALSE)
sweep.stats_brain <- summarizeSweep(sweep.res.list_brain, GT = FALSE)
head(sweep.stats_brain)
sweep.stats_brain[order(sweep.stats_brain$BCreal),]
bcmvn_brain <- find.pK(sweep.stats_brain)
mpK<-as.numeric(as.vector(bcmvn_brain$pK[which.max(bcmvn_brain$BCmetric)]))
DoubletRate = ncol(Seurat_object_P2)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
#DoubletRate = 0.075    # 直接查表，10000细胞对应的doublets rate是～7.6%
#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
Seurat_object_P2 <- FindNeighbors(Seurat_object_P2, reduction = "pca", dims = 1:20)
Seurat_object_P2 <- FindClusters(Seurat_object_P2, resolution = 0.4)
levels(Seurat_object_P2)
#[1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"
Seurat_object_P2 <- RunUMAP(Seurat_object_P2, reduction = "pca", dims = 1:20)
Seurat_object_P2 <- RunTSNE(Seurat_object_P2, reduction = "pca", dims = 1:20)
annotations <- Seurat_object_P2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(DoubletRate*nrow(Seurat_object_P2@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
Seurat_object_P2 <- doubletFinder_v3(Seurat_object_P2, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
Seurat_object_P2 <- doubletFinder_v3(Seurat_object_P2, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.12_14", sct = FALSE)
pdf("DoubletFinder_P2.pdf")
DimPlot(Seurat_object_P2, reduction = "tsne", group.by = "DF.classifications_0.25_0.12_14", pt.size = 0.6)
dev.off()
Seurat_object_P2@meta.data$singledouble<-Seurat_object_P2@meta.data$'DF.classifications_0.25_0.12_14'
Seurat_object_P2.singlet <- subset(Seurat_object_P2, subset = singledouble == "Singlet")
Seurat_object_P2$'DF.classifications_0.25_0.12_14'<-Idents(Seurat_object_P2)
DimPlot(Seurat_object_P2, reduction = "tsne", label = TRUE, pt.size=0.6,label.size = 8, group.by = "DF.classifications_0.25_0.12_14")
Seurat_object_P2_singlet<-Seurat_object_P2.singlet
plot(x=Seurat_object_P2_singlet@meta.data$nCount_RNA,y=Seurat_object_P2_singlet@meta.data$nFeature_RNA)
# check the metadata in the new Seurat objects
head(Seurat_object_P2_singlet@meta.data)
tail(Seurat_object_P2_singlet@meta.data)
# Create .RData object to load at any time
save(Seurat_object_P2_singlet, file="Seurat_object_P2_singlet.RData")
Seurat_object_P2_singlet$log10GenesPerUMI <- log10(Seurat_object_P2_singlet$nFeature_RNA) / log10(Seurat_object_P2_singlet$nCount_RNA)
Seurat_object_P2_singlet$mitoRatio <- PercentageFeatureSet(object = Seurat_object_P2_singlet, pattern = "^MT-")
Seurat_object_P2_singlet$mitoRatio <- Seurat_object_P2_singlet@meta.data$mitoRatio / 100
Seurat_object_P2_singletmetadata <- Seurat_object_P2_singlet@meta.data
Seurat_object_P2_singletmetadata$cells <- rownames(Seurat_object_P2_singletmetadata)
Seurat_object_P2_singletmetadata <- Seurat_object_P2_singletmetadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
Seurat_object_P2_singlet
Seurat_object_P2_singlet@meta.data <- Seurat_object_P2_singletmetadata
counts <- GetAssayData(object = Seurat_object_P2_singlet, slot = "counts")
Seurat_object_P2_singlet <- CreateSeuratObject(counts, meta.data = Seurat_object_P2_singlet@meta.data)
Seurat_object_P2_singlet$label <- "P2"
P2_norm <- NormalizeData(Seurat_object_P2_singlet, normalization.method = "LogNormalize", scale.factor = 10000)
P2_norm <- FindVariableFeatures(P2_norm, selection.method = "vst", nfeatures = 2000)
Seurat_object_P2_singlet$log10GenesPerUMI <- log10(Seurat_object_P2_singlet$nFeature_RNA) / log10(Seurat_object_P2_singlet$nCount_RNA)
Seurat_object_P2_singlet$mitoRatio <- PercentageFeatureSet(object = Seurat_object_P2_singlet, pattern = "^MT-")
Seurat_object_P2_singlet$mitoRatio <- Seurat_object_P2_singlet@meta.data$mitoRatio / 100
VlnPlot(Seurat_object_P2_singlet, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3)
plot1 <- FeatureScatter(Seurat_object_P2_singlet, feature1 = "nCount_RNA", feature2 = "mitoRatio")
plot2 <- FeatureScatter(Seurat_object_P2_singlet, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
filtered_P2 <- subset(x = Seurat_object_P2_singlet, 
                      subset= (nCount_RNA < 8000) & 
                        (nFeature_RNA > 100) &
                        (nFeature_RNA < 4000) & 
                        (log10GenesPerUMI > 0.70) & 
                        (mitoRatio < 0.3))
Seurat_object_P2
#An object of class Seurat 
#9463 features across 1300 samples within 2 assays 
#Active assay: SCT (9458 features, 3000 variable features)
#1 other assay present: RNA
#3 dimensional reductions calculated: pca, umap, tsne
Seurat_object_P2_singlet
#An object of class Seurat 
#9458 features across 1286 samples within 1 assay 
#Active assay: RNA (9458 features, 0 variable features)
filtered_P2
#An object of class Seurat 
#9458 features across 1250 samples within 1 assay 
#Active assay: RNA (9458 features, 0 variable features)
counts <- GetAssayData(object = filtered_P2, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 20
filtered_counts <- counts[keep_genes, ]
filtered_P2 <- CreateSeuratObject(filtered_counts, meta.data = Seurat_object_P2_singlet@meta.data)
filtered_P2$label <- "P2"
filtered_P2
#An object of class Seurat 
#9317 features across 1250 samples within 1 assay 
#Active assay: RNA (9317 features, 0 variable features)
save(filtered_P2, file="filtered_P2.RData")
filtered_P2_norm<-NormalizeData(filtered_P2)

##################第3个病人：GSM3304011_P3_Tumor#########################################
ScRNA_exp_P3 <- read.table(
  "D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/scRNA-seq/GSE117570_RAW/GSM3304011_P3_Tumor_processed_data.txt.gz",
  row.names = 1,
  header = T)
Seurat_object_P3 <- CreateSeuratObject(
  counts = ScRNA_exp_P3, 
  min.cells = 3, 
  min.features = 200)
Seurat_object_P3 <- NormalizeData(Seurat_object_P3)
Seurat_object_P3 <- FindVariableFeatures(Seurat_object_P3, nfeatures = 2000)
Seurat_object_P3[["percent.mt"]] <- PercentageFeatureSet(Seurat_object_P3, pattern = "^MT-")
Seurat_object_P3
Seurat_object_P3 <- ScaleData(Seurat_object_P3)
Seurat_object_P3 <- RunPCA(Seurat_object_P3)
Seurat_object_P3 <- RunUMAP(Seurat_object_P3, reduction = "pca", dims = 1:20)
Seurat_object_P3 <- RunTSNE(Seurat_object_P3, reduction = "pca", dims = 1:20)
## Pre-process Seurat object (sctransform) -----------------------------------------------------------------------------------
Seurat_object_P3 <- SCTransform(Seurat_object_P3)
Seurat_object_P3 <- RunPCA(Seurat_object_P3)
Seurat_object_P3 <- RunUMAP(Seurat_object_P3, reduction = "pca", dims=1:20)
Seurat_object_P3 <- RunTSNE(Seurat_object_P3, reduction = "pca", dims=1:20)
sweep.res.list_brain <- paramSweep_v3(Seurat_object_P3, PCs = 1:20, sct = FALSE)
sweep.stats_brain <- summarizeSweep(sweep.res.list_brain, GT = FALSE)
head(sweep.stats_brain)
sweep.stats_brain[order(sweep.stats_brain$BCreal),]
bcmvn_brain <- find.pK(sweep.stats_brain)
mpK<-as.numeric(as.vector(bcmvn_brain$pK[which.max(bcmvn_brain$BCmetric)]))
DoubletRate = ncol(Seurat_object_P3)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
#DoubletRate = 0.075    # 直接查表，10000细胞对应的doublets rate是～7.6%
#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
Seurat_object_P3 <- FindNeighbors(Seurat_object_P3, reduction = "pca", dims = 1:20)
Seurat_object_P3 <- FindClusters(Seurat_object_P3, resolution = 0.4)
levels(Seurat_object_P3)
#[1] "0"  "1"  "2"
Seurat_object_P3 <- RunUMAP(Seurat_object_P3, reduction = "pca", dims = 1:20)
Seurat_object_P3 <- RunTSNE(Seurat_object_P3, reduction = "pca", dims = 1:20)
annotations <- Seurat_object_P3@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(DoubletRate*nrow(Seurat_object_P3@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
Seurat_object_P3 <- doubletFinder_v3(Seurat_object_P3, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
Seurat_object_P3 <- doubletFinder_v3(Seurat_object_P3, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.05_1", sct = FALSE)
pdf("DoubletFinder_P3.pdf")
DimPlot(Seurat_object_P3, reduction = "tsne", group.by = "DF.classifications_0.25_0.05_1", pt.size = 0.6)
dev.off()
Seurat_object_P3@meta.data$singledouble<-Seurat_object_P3@meta.data$'DF.classifications_0.25_0.05_1'
Seurat_object_P3.singlet <- subset(Seurat_object_P3, subset = singledouble == "Singlet")
Seurat_object_P3$'DF.classifications_0.25_0.05_1'<-Idents(Seurat_object_P3)
DimPlot(Seurat_object_P3, reduction = "tsne", label = TRUE, pt.size=0.6,label.size = 8, group.by = "DF.classifications_0.25_0.05_1")
Seurat_object_P3_singlet<-Seurat_object_P3.singlet
plot(x=Seurat_object_P3_singlet@meta.data$nCount_RNA,y=Seurat_object_P3_singlet@meta.data$nFeature_RNA)
# check the metadata in the new Seurat objects
head(Seurat_object_P3_singlet@meta.data)
tail(Seurat_object_P3_singlet@meta.data)
# Create .RData object to load at any time
save(Seurat_object_P3_singlet, file="Seurat_object_P3_singlet.RData")
Seurat_object_P3_singlet$log10GenesPerUMI <- log10(Seurat_object_P3_singlet$nFeature_RNA) / log10(Seurat_object_P3_singlet$nCount_RNA)
Seurat_object_P3_singlet$mitoRatio <- PercentageFeatureSet(object = Seurat_object_P3_singlet, pattern = "^MT-")
Seurat_object_P3_singlet$mitoRatio <- Seurat_object_P3_singlet@meta.data$mitoRatio / 100
Seurat_object_P3_singletmetadata <- Seurat_object_P3_singlet@meta.data
Seurat_object_P3_singletmetadata$cells <- rownames(Seurat_object_P3_singletmetadata)
Seurat_object_P3_singletmetadata <- Seurat_object_P3_singletmetadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
Seurat_object_P3_singlet
Seurat_object_P3_singlet@meta.data <- Seurat_object_P3_singletmetadata
counts <- GetAssayData(object = Seurat_object_P3_singlet, slot = "counts")
Seurat_object_P3_singlet <- CreateSeuratObject(counts, meta.data = Seurat_object_P3_singlet@meta.data)
Seurat_object_P3_singlet$label <- "P3"
P3_norm <- NormalizeData(Seurat_object_P3_singlet, normalization.method = "LogNormalize", scale.factor = 10000)
P3_norm <- FindVariableFeatures(P3_norm, selection.method = "vst", nfeatures = 2000)
Seurat_object_P3_singlet$log10GenesPerUMI <- log10(Seurat_object_P3_singlet$nFeature_RNA) / log10(Seurat_object_P3_singlet$nCount_RNA)
Seurat_object_P3_singlet$mitoRatio <- PercentageFeatureSet(object = Seurat_object_P3_singlet, pattern = "^MT-")
Seurat_object_P3_singlet$mitoRatio <- Seurat_object_P3_singlet@meta.data$mitoRatio / 100
VlnPlot(Seurat_object_P3_singlet, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3)
plot1 <- FeatureScatter(Seurat_object_P3_singlet, feature1 = "nCount_RNA", feature2 = "mitoRatio")
plot2 <- FeatureScatter(Seurat_object_P3_singlet, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
filtered_P3 <- subset(x = Seurat_object_P3_singlet, 
                      subset= (nCount_RNA < 12600) & 
                        (nFeature_RNA > 100) &
                        (nFeature_RNA < 2500) & 
                        (log10GenesPerUMI > 0.70) & 
                        (mitoRatio < 0.3))
Seurat_object_P3
#An object of class Seurat 
#3211 features across 328 samples within 2 assays 
#Active assay: SCT (3211 features, 3000 variable features)
#1 other assay present: RNA
#3 dimensional reductions calculated: pca, umap, tsne
Seurat_object_P3_singlet
#An object of class Seurat 
#3211 features across 327 samples within 1 assay 
#Active assay: RNA (3211 features, 0 variable features)
filtered_P3
#An object of class Seurat 
#3211 features across 319 samples within 1 assay 
#Active assay: RNA (3211 features, 0 variable features)
counts <- GetAssayData(object = filtered_P3, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 20
filtered_counts <- counts[keep_genes, ]
filtered_P3 <- CreateSeuratObject(filtered_counts, meta.data = Seurat_object_P3_singlet@meta.data)
filtered_P3$label <- "P3"
filtered_P3
#An object of class Seurat 
#3205 features across 319 samples within 1 assay 
#Active assay: RNA (3205 features, 0 variable features)
save(filtered_P3, file="filtered_P3.RData")
filtered_P3_norm<-NormalizeData(filtered_P3)

##################第4个病人：GSM3304013_P4_Tumor#########################################
ScRNA_exp_P4 <- read.table(
  "D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/scRNA-seq/GSE117570_RAW/GSM3304013_P4_Tumor_processed_data.txt.gz",
  row.names = 1,
  header = T)
Seurat_object_P4 <- CreateSeuratObject(
  counts = ScRNA_exp_P4, 
  min.cells = 3, 
  min.features = 200)
Seurat_object_P4 <- NormalizeData(Seurat_object_P4)
Seurat_object_P4 <- FindVariableFeatures(Seurat_object_P4, nfeatures = 2000)
Seurat_object_P4[["percent.mt"]] <- PercentageFeatureSet(Seurat_object_P4, pattern = "^MT-")
Seurat_object_P4
Seurat_object_P4 <- ScaleData(Seurat_object_P4)
Seurat_object_P4 <- RunPCA(Seurat_object_P4)
Seurat_object_P4 <- RunUMAP(Seurat_object_P4, reduction = "pca", dims = 1:20)
Seurat_object_P4 <- RunTSNE(Seurat_object_P4, reduction = "pca", dims = 1:20)
## Pre-process Seurat object (sctransform) -----------------------------------------------------------------------------------
Seurat_object_P4 <- SCTransform(Seurat_object_P4)
Seurat_object_P4 <- RunPCA(Seurat_object_P4)
Seurat_object_P4 <- RunUMAP(Seurat_object_P4, reduction = "pca", dims=1:20)
Seurat_object_P4 <- RunTSNE(Seurat_object_P4, reduction = "pca", dims=1:20)
sweep.res.list_brain <- paramSweep_v3(Seurat_object_P4, PCs = 1:20, sct = FALSE)
sweep.stats_brain <- summarizeSweep(sweep.res.list_brain, GT = FALSE)
head(sweep.stats_brain)
sweep.stats_brain[order(sweep.stats_brain$BCreal),]
bcmvn_brain <- find.pK(sweep.stats_brain)
mpK<-as.numeric(as.vector(bcmvn_brain$pK[which.max(bcmvn_brain$BCmetric)]))
DoubletRate = ncol(Seurat_object_P4)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
#DoubletRate = 0.075    # 直接查表，10000细胞对应的doublets rate是～7.6%
#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
Seurat_object_P4 <- FindNeighbors(Seurat_object_P4, reduction = "pca", dims = 1:20)
Seurat_object_P4 <- FindClusters(Seurat_object_P4, resolution = 0.4)
levels(Seurat_object_P4)
#[1] "0"  "1"  "2"  "3"  "4"
Seurat_object_P4 <- RunUMAP(Seurat_object_P4, reduction = "pca", dims = 1:20)
Seurat_object_P4 <- RunTSNE(Seurat_object_P4, reduction = "pca", dims = 1:20)
annotations <- Seurat_object_P4@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(DoubletRate*nrow(Seurat_object_P4@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
Seurat_object_P4 <- doubletFinder_v3(Seurat_object_P4, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
Seurat_object_P4 <- doubletFinder_v3(Seurat_object_P4, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.29_16", sct = FALSE)
pdf("DoubletFinder_P4.pdf")
DimPlot(Seurat_object_P4, reduction = "tsne", group.by = "DF.classifications_0.25_0.29_16", pt.size = 0.6)
dev.off()
Seurat_object_P4@meta.data$singledouble<-Seurat_object_P4@meta.data$'DF.classifications_0.25_0.29_16'
Seurat_object_P4.singlet <- subset(Seurat_object_P4, subset = singledouble == "Singlet")
Seurat_object_P4$'DF.classifications_0.25_0.29_16'<-Idents(Seurat_object_P4)
DimPlot(Seurat_object_P4, reduction = "tsne", label = TRUE, pt.size=0.6,label.size = 8, group.by = "DF.classifications_0.25_0.29_16")
Seurat_object_P4_singlet<-Seurat_object_P4.singlet
plot(x=Seurat_object_P4_singlet@meta.data$nCount_RNA,y=Seurat_object_P4_singlet@meta.data$nFeature_RNA)
# check the metadata in the new Seurat objects
head(Seurat_object_P4_singlet@meta.data)
tail(Seurat_object_P4_singlet@meta.data)
# Create .RData object to load at any time
save(Seurat_object_P4_singlet, file="Seurat_object_P4_singlet.RData")
Seurat_object_P4_singlet$log10GenesPerUMI <- log10(Seurat_object_P4_singlet$nFeature_RNA) / log10(Seurat_object_P4_singlet$nCount_RNA)
Seurat_object_P4_singlet$mitoRatio <- PercentageFeatureSet(object = Seurat_object_P4_singlet, pattern = "^MT-")
Seurat_object_P4_singlet$mitoRatio <- Seurat_object_P4_singlet@meta.data$mitoRatio / 100
Seurat_object_P4_singletmetadata <- Seurat_object_P4_singlet@meta.data
Seurat_object_P4_singletmetadata$cells <- rownames(Seurat_object_P4_singletmetadata)
Seurat_object_P4_singletmetadata <- Seurat_object_P4_singletmetadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
Seurat_object_P4_singlet
Seurat_object_P4_singlet@meta.data <- Seurat_object_P4_singletmetadata
counts <- GetAssayData(object = Seurat_object_P4_singlet, slot = "counts")
Seurat_object_P4_singlet <- CreateSeuratObject(counts, meta.data = Seurat_object_P4_singlet@meta.data)
Seurat_object_P4_singlet$label <- "P4"
P4_norm <- NormalizeData(Seurat_object_P4_singlet, normalization.method = "LogNormalize", scale.factor = 10000)
P4_norm <- FindVariableFeatures(P4_norm, selection.method = "vst", nfeatures = 2000)
Seurat_object_P4_singlet$log10GenesPerUMI <- log10(Seurat_object_P4_singlet$nFeature_RNA) / log10(Seurat_object_P4_singlet$nCount_RNA)
Seurat_object_P4_singlet$mitoRatio <- PercentageFeatureSet(object = Seurat_object_P4_singlet, pattern = "^MT-")
Seurat_object_P4_singlet$mitoRatio <- Seurat_object_P4_singlet@meta.data$mitoRatio / 100
VlnPlot(Seurat_object_P4_singlet, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3)
plot1 <- FeatureScatter(Seurat_object_P4_singlet, feature1 = "nCount_RNA", feature2 = "mitoRatio")
plot2 <- FeatureScatter(Seurat_object_P4_singlet, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
filtered_P4 <- subset(x = Seurat_object_P4_singlet, 
                      subset= (nCount_RNA < 9000) & 
                        (nFeature_RNA > 100) &
                        (nFeature_RNA < 3500) & 
                        (log10GenesPerUMI > 0.70) & 
                        (mitoRatio < 0.3))
Seurat_object_P4
#An object of class Seurat 
#7492 features across 1423 samples within 2 assays 
#Active assay: SCT (7486 features, 3000 variable features)
#1 other assay present: RNA
#3 dimensional reductions calculated: pca, umap, tsne
Seurat_object_P4_singlet
#An object of class Seurat 
#7486 features across 1407 samples within 1 assay 
#Active assay: RNA (7486 features, 0 variable features)
filtered_P4
#An object of class Seurat 
#7486 features across 1364 samples within 1 assay 
#Active assay: RNA (7486 features, 0 variable features)
counts <- GetAssayData(object = filtered_P4, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 20
filtered_counts <- counts[keep_genes, ]
filtered_P4 <- CreateSeuratObject(filtered_counts, meta.data = Seurat_object_P4_singlet@meta.data)
filtered_P4$label <- "P4"
filtered_P4
#An object of class Seurat 
#7454 features across 1364 samples within 1 assay 
#Active assay: RNA (7454 features, 0 variable features)
save(filtered_P4, file="filtered_P4.RData")
filtered_P4_norm<-NormalizeData(filtered_P4)

###############################################merged samples####################################################################
NSCLC.anchors <- FindIntegrationAnchors(object.list = list(filtered_P1_norm, filtered_P2_norm, filtered_P3_norm, 
                                                           filtered_P4_norm), dims = 1:30)
save(NSCLC.anchors, file="integrated.anchors_seurat_4samples.RData")
NSCLC.combined <- IntegrateData(anchorset = NSCLC.anchors, dims = 1:30)
NSCLC.combined <- FindVariableFeatures(NSCLC.combined, selection.method = "vst", nfeatures = 2000)
save(NSCLC.combined, file="integrated.combined_seurat_4samples.RData")

pdf("VlnPlot_NSCLC_4patients.pdf",height = 4, width = 6)
VlnPlot(NSCLC.combined, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), pt.size = 0, group.by = "label", ncol = 3)
dev.off()
pdf("Correlation_NSCLC_4patients.pdf",height = 4, width = 6)
FeatureScatter(NSCLC.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "label")
dev.off()
DefaultAssay(NSCLC.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
NSCLC.combined <- ScaleData(NSCLC.combined, verbose = FALSE, vars.to.regress = c("nUMI", "mitoRatio"))
#Scaling is an essential step in the Seurat workflow, but only on genes that will be used as input to PCA. Therefore, the default in ScaleData() is only to perform scaling on the previously identified variable features (2,000 by default). 
NSCLC.combined <- ScaleData(NSCLC.combined)
save(NSCLC.combined, file="NSCLC.combined_scaled.RData")
NSCLC.combined <- RunPCA(NSCLC.combined, npcs = 30, verbose = FALSE)
print(NSCLC.combined[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(NSCLC.combined, dims = 1:2, reduction = "pca")
DimPlot(NSCLC.combined, reduction = "pca")
DimHeatmap(NSCLC.combined, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(NSCLC.combined, dims = 1:30, cells = 500, balanced = TRUE)
NSCLC.combined <- JackStraw(NSCLC.combined, num.replicate = 100, dims = 30)
NSCLC.combined <- ScoreJackStraw(NSCLC.combined, dims = 1:30)
library("cowplot")
pdf("Determine_dimensionality.pdf", width = 24, height = 18)
p1 <- JackStrawPlot(NSCLC.combined, dims = 1:30)
p2 <- ElbowPlot(NSCLC.combined,ndims = 30)
plot_grid(p1, p2)
dev.off()

################################### determine the number of pcs##############################################
#筛选标准
#1.主成分累积贡献大于90%
#2.PC本身对方差贡献小于5%
#3.两个连续PCs之间差异小于0.1%
#将要检测的seurat对象传递给sce
library("ggplot2")
sce=NSCLC.combined
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
#[1] 14
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
save(NSCLC.combined, file="integrated.combined_beforepcs.RData")


################################ determine the resolution###########################################################
library(Seurat)
library(clustree)
library(SeuratData)
DefaultAssay(NSCLC.combined) <- "RNA"
NSCLC.combined <- FindVariableFeatures(NSCLC.combined, selection.method = "vst", nfeatures = 2000)
NSCLC.combined <- ScaleData(NSCLC.combined , verbose = FALSE)
NSCLC.combined <- RunPCA(NSCLC.combined)
NSCLC.combined <- RunUMAP(NSCLC.combined, reduction = "pca", dims = 1:25)
NSCLC.combined <- RunTSNE(NSCLC.combined, reduction = "pca", dims = 1:25)
NSCLC.combined <- FindNeighbors(NSCLC.combined, reduction = "pca", dims = 1:25)
NSCLC.combined <- FindClusters(
  object = NSCLC.combined,
  resolution = c(seq(0,1.6,.2))
)
clustree(NSCLC.combined@meta.data, prefix = "RNA_snn_res.")
pdf("Determine_resolution.pdf")
clustree(NSCLC.combined, prefix = "RNA_snn_res.") +
  guides(edge_colour = FALSE, edge_alpha = FALSE) +scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "blue", high = "red")+
  theme(legend.position = "bottom")
dev.off()
r=0.6
NSCLC.combined <- FindClusters(NSCLC.combined, resolution = r)
levels(NSCLC.combined)
#[1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20"
save(NSCLC.combined,file="NSCLC.combined.res0.6.RData")
NSCLC.pca14.markers <- FindAllMarkers(object = NSCLC.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                                      test.use = "wilcox")
save(NSCLC.pca14.markers,file="NSCLC.pca14.markers.res0.6.RData")
library(dplyr)
top30<-NSCLC.pca14.markers %>% group_by(cluster) %>% top_n(n=30,wt=avg_log2FC)
write.table(top30,"top30.pca14.res0.6.markers.csv",sep=",",quote=F)
head(Idents(NSCLC.combined), 5)
write.table(NSCLC.pca14.markers,"NSCLC.pca14.markers.csv",sep=",",quote=F)


################################细胞亚群注释#####################################################################
library("SingleR")
library("celldex")
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se
DefaultAssay(NSCLC.combined)<-"integrated"
testdata <- GetAssayData(NSCLC.combined, slot="data")
clusters <- NSCLC.combined@meta.data$seurat_clusters
cellpred <- SingleR(test=testdata, ref=hpca.se, labels=hpca.se$label.main, method = "cluster", clusters = clusters, 
                    de.method="wilcox")
table(cellpred$labels)
#B_cell Endothelial_cells  Epithelial_cells       Fibroblasts        Macrophage          Monocyte           NK_cell 
#     1                 2                 7                 1                 4                 3                 3 
pdf("SingleR_test.pdf")
plotScoreHeatmap(cellpred)
dev.off()
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"celltype_singleR.csv",row.names = F)
NSCLC.combined@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  NSCLC.combined@meta.data[which(NSCLC.combined@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

p1 = DimPlot(NSCLC.combined, group.by="celltype", label=T, label.size=5, reduction='tsne')
ggsave("NSCLC.combined.tSNE_celltype.pdf", p1, width=8 ,height=6)
# annotation
new.cluster.ids<-celltype$celltype
names(new.cluster.ids) <- levels(NSCLC.combined)
NSCLC.combined <- RenameIdents(NSCLC.combined, new.cluster.ids)
NSCLC.combined$celltype<-Idents(NSCLC.combined)
save(NSCLC.combined,file="NSCLC.combined.pca14.res0.6.afteranno.RData")
NSCLC.combined.markers <- FindAllMarkers(object = NSCLC.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                                         test.use = "wilcox")
write.table(NSCLC.combined.markers,"NSCLC.combined.markers_afterAnnotation.csv",sep=",",quote=F)
top30<-NSCLC.combined.markers %>% group_by(cluster) %>% top_n(n=30,wt=avg_log2FC)
write.table(top30,"top30.NSCLC.combined.markers.csv",sep=",",quote=F)
save(NSCLC.combined.markers,file="NSCLC.combined.markers.afteranno.RData")
# calculate the percentage and count the numbers
prop.table(table(Idents(NSCLC.combined), NSCLC.combined$label))
allsampleprop.each <-as.data.frame(prop.table(table(Idents(NSCLC.combined), NSCLC.combined$label)))
prop.table(table(Idents(NSCLC.combined)))
allsampleprop.total <-as.data.frame(prop.table(table(Idents(NSCLC.combined))))
write.csv(x = allsampleprop.each,file = 'anno.allsample.each.prop.csv',quote = T,row.names = T)
write.csv(x = allsampleprop.total,file = 'anno.allsample.total.prop.csv',quote = T,row.names = T)
table(Idents(NSCLC.combined))
#Epithelial_cells          Monocyte           NK_cell        Macrophage            B_cell       Fibroblasts Endothelial_cells 
#            1815               957               676               980               170                52                71 
pro.total <- table(Idents(NSCLC.combined),NSCLC.combined$label)
table(Idents(NSCLC.combined),NSCLC.combined$label)
pro.each <- table(Idents(NSCLC.combined),NSCLC.combined$label)
write.csv(x =pro.total,file = 'anno.pro.total.csv',quote = T,row.names = T)
write.csv(x =pro.each,file = 'anno.pro.each.csv',quote = T,row.names = T)
# visualization
NSCLC.combined$celltype <- Idents(NSCLC.combined)
pdf(paste0("tsne.pca14.res",r,".splitbyLabel_NoLegend.pdf"),width=40,height=10)
DimPlot(NSCLC.combined, reduction = "tsne", label = TRUE, pt.size=0.6, label.size = 0, split.by = 'label', group.by = "celltype")+
  NoLegend()
dev.off()
pdf(paste0("tsne.pca14.res",r,".splitbyLabel_Legend.pdf"),width=40,height=10)
DimPlot(NSCLC.combined, reduction = "tsne", label = TRUE, pt.size=0.6,label.size = 8, split.by = 'label', group.by = 'celltype')
dev.off()
pdf(paste0("tsne.pca14.res",r,"_Legend.pdf"),width=10,height=10)
DimPlot(NSCLC.combined, reduction = "tsne", label = TRUE, pt.size=0.6,label.size = 8, group.by = 'celltype')
dev.off()
pdf(paste0("tsne.pca14.res",r,"_NoLegend.pdf"),width=10,height=10)
DimPlot(NSCLC.combined, reduction = "tsne", label = FALSE, pt.size=0.6,label.size = 8, group.by = 'celltype')+NoLegend()
dev.off()
pdf("tsne.merged.pdf",width=10,height=10)
DimPlot(NSCLC.combined, reduction = "tsne", label = TRUE, pt.size=0.6,label.size = 8, group.by = "label")+NoLegend()
dev.off()

###########################################################################根据已有的定义判断细胞是否为免疫细胞####################################################################
load("NSCLC.combined.pca14.res0.6.afteranno.RData")
levels(NSCLC.combined)
#[1] "Epithelial_cells"  "Monocyte"          "NK_cell"           "Macrophage"        "B_cell"
#     "Fibroblasts"       "Endothelial_cells"
r=0.6
pdf(paste0("tsne.pca14.res",r,".integrated.celltypes.pdf"),width=10,height=10)
DimPlot(NSCLC.combined, reduction = "tsne", label = TRUE, pt.size=0.6,label.size = 8, group.by = "celltype")+NoLegend()
dev.off()
immune.labels<-c("Non-immune", "Immune", "Immune", "Immune", "Immune", "Non-immune", "Non-immune")
table(immune.labels)
#immune.labels
#Immune Non-immune 
#4          3 
names(immune.labels) <- levels(NSCLC.combined)
NSCLC.combined <- RenameIdents(NSCLC.combined, immune.labels)
levels(NSCLC.combined)
#[1] "Non-immune" "Immune"
NSCLC.combined$immune<-Idents(NSCLC.combined)
head(NSCLC.combined@meta.data)
save(NSCLC.combined,file="NSCLC.combined.immuneORnonimmuneClass.RData")
pdf(paste0("tsne.pca14.res",r,".Immune.pdf"),width=10,height=10)
DimPlot(NSCLC.combined, reduction = "tsne", label = TRUE, pt.size=0.6,label.size = 8, group.by = "immune")+NoLegend()
dev.off()

#########################################热图可视化##################################################################
load("NSCLC.combined.pca14.res0.6.afteranno.RData")
pdf("heatmap.top30.pdf",width=24,height=18)
DoHeatmap(NSCLC.combined,features=top30$gene,cells = 1:500, size = 4, angle = 90, disp.min=-2, disp.max=2) + scale_fill_gradientn(colours=c("blue","white","red"))
dev.off()

#########################################气泡图可视化#################################################################
#从panglaodb数据库中挑选相应亚群的marker基因
features.plot <- c("MUC1","TSTD1","ELF3",
                   "ANPEP","HBEGF","ICAM3",
                   "CD69","IL32","CXCR4",
                   "CXCL16","S100A4","SLC11A1",
                   "CD40","SWAP70","MEF2C",
                   "EGR1","FOSB","CBX1",
                   "TGFBR2","IL6ST","IFI16")
pdf("NSCLC.dittoDotPlot.pdf",width = 11, height = 6)
DotPlot(object = NSCLC.combined, features=features.plot,dot.scale = 6,cols = c("lightgrey", "red")) + RotatedAxis()
dev.off()

#####################################tsne可视化图各个亚群#############################################################
levels(NSCLC.combined)
#[1] "Epithelial_cells"  "Monocyte"          "NK_cell"           "Macrophage"        "B_cell"            "Fibroblasts"      
#[7] "Endothelial_cells"

pdf("subpopulations/Epithelial_cells_tsne.pdf", width = 10, height = 10)
DimPlot(NSCLC.combined, reduction = "tsne",group.by ="celltype",pt.size = 0.8,
        cols=c('Epithelial_cells'='red','Monocyte'='grey','NK_cell'='grey','Macrophage'='grey','B_cell'='grey',
               'Fibroblasts'='grey','Endothelial_cells'='grey'))+ NoLegend()
dev.off()

pdf("subpopulations/Monocyte_tsne.pdf", width = 10, height = 10)
DimPlot(NSCLC.combined, reduction = "tsne",group.by ="celltype",pt.size = 0.8,
        cols=c('Epithelial_cells'='grey','Monocyte'='red','NK_cell'='grey','Macrophage'='grey','B_cell'='grey',
               'Fibroblasts'='grey','Endothelial_cells'='grey'))+ NoLegend()
dev.off()

pdf("subpopulations/NK_cell_tsne.pdf", width = 10, height = 10)
DimPlot(NSCLC.combined, reduction = "tsne",group.by ="celltype",pt.size = 0.8,
        cols=c('Epithelial_cells'='grey','Monocyte'='grey','NK_cell'='red','Macrophage'='grey','B_cell'='grey',
               'Fibroblasts'='grey','Endothelial_cells'='grey'))+ NoLegend()
dev.off()

pdf("subpopulations/Macrophage_tsne.pdf", width = 10, height = 10)
DimPlot(NSCLC.combined, reduction = "tsne",group.by ="celltype",pt.size = 0.8,
        cols=c('Epithelial_cells'='grey','Monocyte'='grey','NK_cell'='grey','Macrophage'='red','B_cell'='grey',
               'Fibroblasts'='grey','Endothelial_cells'='grey'))+ NoLegend()
dev.off()

pdf("subpopulations/B_cell_tsne.pdf", width = 10, height = 10)
DimPlot(NSCLC.combined, reduction = "tsne",group.by ="celltype",pt.size = 0.8,
        cols=c('Epithelial_cells'='grey','Monocyte'='grey','NK_cell'='grey','Macrophage'='grey','B_cell'='red',
               'Fibroblasts'='grey','Endothelial_cells'='grey'))+ NoLegend()
dev.off()

pdf("subpopulations/Fibroblasts_tsne.pdf", width = 10, height = 10)
DimPlot(NSCLC.combined, reduction = "tsne",group.by ="celltype",pt.size = 0.8,
        cols=c('Epithelial_cells'='grey','Monocyte'='grey','NK_cell'='grey','Macrophage'='grey','B_cell'='grey',
               'Fibroblasts'='red','Endothelial_cells'='grey'))+ NoLegend()
dev.off()

pdf("subpopulations/Endothelial_cells_tsne.pdf", width = 10, height = 10)
DimPlot(NSCLC.combined, reduction = "tsne",group.by ="celltype",pt.size = 0.8,
        cols=c('Epithelial_cells'='grey','Monocyte'='grey','NK_cell'='grey','Macrophage'='grey','B_cell'='grey',
               'Fibroblasts'='grey','Endothelial_cells'='red'))+ NoLegend()
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
ids=bitr(NSCLC.combined.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
sce.markers=merge(NSCLC.combined.markers,ids,by.x='gene',by.y='SYMBOL')
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
setwd("D:/MyProjects/scRNA_immune/scRNA_CAF_NSCLC/01_scRNA/subpopulations/EnrichmentAnalysis/")
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

