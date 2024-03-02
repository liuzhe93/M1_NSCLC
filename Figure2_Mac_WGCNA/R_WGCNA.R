setwd("/Users/liuzhe/Desktop/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/WGCNA")
remove(list=ls())

#ref: https://cloud.tencent.com/developer/article/1498123

#clinical data processing
rt<-read.csv("/Users/liuzhe/Desktop/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv")
var="Macrophages.M1"
rt=rt[,c("sample_name","futime","fustat",var)]
group=ifelse(rt[,4]>median(rt[,4]),"High","Low")
rt$group<-group


#expression profile data processing
exp<-read.csv("/Users/liuzhe/Desktop/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ExprData/merged_expr_dat_remove_batch_effect.csv",header = T, row.names = 1)
exp_t<-t(exp)
exp_t<-as.data.frame(exp_t)
exp_t$sample_name<-row.names(exp_t)
exp_t$sample_name<-substr(exp_t$sample_name, 1, 12)
exp_t$sample_name<-gsub("\\.","-",exp_t$sample_name)

#expression profile with group label
#################################step1_构建WGCNA所需表达矩阵#####################################################################
merged_data<-merge(exp_t, rt, by = "sample_name")
val1=c("futime", "fustat", "Macrophages.M1")
tmp_data<-merged_data[, !colnames(merged_data) %in% val1]
group_list<-factor(tmp_data$group)
table(group_list)
#group_list
#High  Low 
#806  967
for(i in 1:1773){
  rownames(tmp_data)[i]<-paste0(tmp_data$sample_name,"_",i)
}
tmp_data<-tmp_data[,-1]
tmp_data<-tmp_data[,-ncol(tmp_data)]
datTraits <- data.frame(row.names=rownames(tmp_data),
                        subtype=group_list)
library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
datExpr <- tmp_data
nSamples = nrow(datExpr)
nSamples
#[1] 1773
datExpr_t<-t(datExpr)

###########基因水平的过滤
## 筛选中位绝对偏差前75%的基因，至少MAD大于0.01
## 筛选后会降低运算量，也会失去部分信息
## 也可不做筛选，使MAD大于0即可
m.mad <- apply(datExpr_t,1,mad)
dataExprVar <- datExpr_t[which(m.mad > 
                                max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
## 转换为样品在行，基因在列的矩阵
dataExpr <- as.data.frame(t(dataExprVar))
dim(dataExpr)
#[1]  1773 11697
dim(datExpr_t)
#[1] 15596  1773
## 检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3)

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)
#[1]  1773 11697
head(dataExpr)[,1:8]
## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")
par(mfrow = c(1,1))
pdf("SampleTree.pdf",width = 20,height = 5)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()


#################################step2_找到合适的beta值#####################################################################
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)
RpowerTable=pickSoftThreshold(datExpr, powerVector=powers)[[2]]
cex1=0.7
pdf(file="softThresholding.pdf")
par(mfrow=c(1,2))
plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2], labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
plot(RpowerTable[,1], RpowerTable[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n")
text(RpowerTable[,1], RpowerTable[,5], labels=powers, cex=cex1,col="red")
dev.off()
power = sft$powerEstimate
power
#[1] 3

##################################step3_对基因聚类成模块并且可视化#############################################################
##一步法网络构建：One-step network construction and module detection##
# power: 上一步计算的软阈值
# maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)；
#  4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可
#  以处理3万个
#  计算资源允许的情况下最好放在一个block里面。
# corType: pearson or bicor
# numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
# saveTOMs：最耗费时间的计算，存储起来，供后续使用
# mergeCutHeight: 合并模块的阈值，越大模块越少
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 100,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0("NSCLC", ".tom"),
                       verbose = 3)
# 根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
# **0 (grey)**表示**未**分入任何模块的基因。 
table(net$colors)
#   0    1    2    3    4    5    6    7    8    9   10   11 
#2950 2880 1363 1179  967  873  532  234  222  197  159  141

## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
table(moduleColors)
#moduleColors
#black        blue       brown       green greenyellow        grey     magenta        pink      purple         red   turquoise      yellow 
#  234        1363        1179         873         141        2950         197         222         159         532        2880         967 
# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
pdf("cluster_dendrograms.pdf")
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


#######################step4_将样本信息添加进去,分析样本性状与模块的关系########################################
design=model.matrix(~0+ datTraits$subtype)
design = as.data.frame(design)
colnames(design)=levels(datTraits$subtype)
moduleColors <- labels2colors(moduleColors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(dataExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, design , use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot
pdf("moduleTraitCor.pdf")
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


############################step5_感兴趣性状的模块的具体基因分析###############################################
## 首先计算模块与基因的相关性矩阵
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="")

## 再计算性状与基因的相关性矩阵
##1-->low; 0-->high
M1 = as.data.frame(design[,2])
names(M1) = "M1"
geneTraitSignificance = as.data.frame(cor(dataExpr, M1, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(M1), sep="")
names(GSPvalue) = paste("p.GS.", names(M1), sep="")

## 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
module = "pink"
module
## [1] "pink"
column = match(module, modNames);
moduleGenes = moduleColors==module;


##############################step6_网络的可视化#####################################################
# -----------------------------------------------------------------------------------
## 首先针对所有基因画热图
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
geneTree = net$dendrograms[[1]]; 
dissTOM = 1-TOMsimilarityFromExpr(dataExpr, power = 3); 
plotTOM = dissTOM^7; 
diag(plotTOM) = NA; 
# TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

## 最后画模块和性状的关系
# Recalculate module eigengenes
MEs = moduleEigengenes(dataExpr, moduleColors)$eigengenes

## 这里把是否属于 M1 表型这个变量用0,1进行数值化。
M1 = as.data.frame(design[,2])
names(M1) = "M1"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, M1))
# Plot the relationships among the eigengenes and the trait

plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), 
                      marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)
# Plot the dendrogram

## 模块与性状的聚类图
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)

## 性状与模块热图
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

################################step7_提取指定模块的基因名################################################
# --------------------------------------------------------------------------------------
# Select module
module = "pink"

# Select module probes
probes = colnames(dataExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
head(modProbes)
## [1] "ACAP1"    "ADAM19"   "ADAMDEC1" "AKNA"     "APOBEC3G" "ARHGAP25"
#######################Step9: 模块的导出##########################################
# --------------------------------------------------------------------------------------
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(dataExpr, power = 3); 
module = "pink";
# Select module probes
probes = colnames(dataExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
## 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
## 模块对应的基因关系矩阵 
cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes, 
  nodeAttr = moduleColors[inModule]
)
save(cyt, file = "Cyt.RData")

