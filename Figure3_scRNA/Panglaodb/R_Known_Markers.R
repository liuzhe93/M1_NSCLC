setwd("D:/MyProjects/scRNA_immune/Mac_Lung/Figure3_scRNA/Panglaodb")
rm(list=ls())

mydata <- read.table("PanglaoDB_markers_27_Mar_2020.tsv", header = T, sep = "\t")
colnames(mydata)<-c("species", "official_gene_symbol", "cell_type", "nicknames", "ubiquitousness_index",
                    "product_description", "gene_type", "canonical_marker", "germ_layer", "organ",
                    "sensitivity_human", "sensitivity_mouse", "specificity_human", "specificity_mouse")

mydata[which(mydata$cell_type == "Monocyte"),]



df[which(df$ID == 4),]


mydata_Epithelial_cells = matrix()
for(i in 1:nrow(mydata)){
  if(mydata$cell_type[i] == "Epithelial cells"){
    tmp = mydata[i,]
    mydata_Epithelial_cells <- rbind(mydata_Epithelial_cells,tmp)
  }
}


mydata_Epithelial_cells <- mydata[mydata$cell_type == "Epithelial cells"]
mydata_Epithelial_cells


cols_name <- c("Epithelial cells", "Monocyte", "NK cell", "Macrophage", "B cell", "Fibroblasts", "Endothelial cells")


"Epithelial_cells"
MUC1
TSTD1
ELF3
CXCL17
AQP3


"Monocyte"
ANPEP
HBEGF


"NK_cell"
CD69
IL32

"Macrophage"
CD14
FCGR3A
CD86
S100A8
CXCL16
CD83
MS4A4A
MS4A6A
RAB20
FMNL1
C5AR1
CD300A
IL1B
MYO1G
SAMSN1
NR4A3
S100A4
MS4A7
SLC11A1
CCL3
CLEC7A
CD74
CSF1R
CD68
FGR
CYBB
ITGAX
TYROBP
RGS1
AIF1
SLC7A7
RUNX3

"B_cell"
CD40
SWAP70


"Fibroblasts"
EGR1
FOSB

"Endothelial_cells"
TGFBR2
