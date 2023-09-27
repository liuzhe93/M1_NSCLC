setwd("D:/MyProjects/scRNA_immune/Mac_Lung/Figure5_M1_prognosis")
rm(list=ls())

library("ConsensusClusterPlus")

lassoSigExp<-read.table("lassoSigExp.txt",header = T, row.names = 1)
d<-lassoSigExp[,3:ncol(lassoSigExp)]
d_t<-t(d)
mads=apply(d_t,1,mad)# mad(x) 绝对中位数差 按行（1）取d数据的中位数
mydata = sweep(d_t,1, apply(d_t,1,median,na.rm=T))
#按行减去中位数，r语言中使用sweep(x, MARGIN, STATS, FUN="-", ...) 对矩阵进行运算。MARGIN为1，表示行的方向上进行运算，
#为2表示列的方向上运算。STATS是运算的参数。FUN为运算函数，默认是减法。
mydata<-as.matrix(mydata)
results <- ConsensusClusterPlus(d = mydata, # 分析矩阵 
                                maxK = 20, # 最大聚类数目 
                                reps = 1000, # 重抽样的次数 
                                pItem = 0.8, # 样品的重抽样比例 
                                clusterAlg = "km", # 使用的聚类算法，可以选择"hc"(hclust), "pam", "km"(k-means) hc对应pearson
                                innerLinkage = "ward.D2", finalLinkage = "ward.D2", 
                                distance = "euclidean", # 计算距离的方法，可以选择pearson、spearman、euclidean、binary/maximum、canberra、minkowski 
                                seed = 123456, # 设置随机种子，方便重复 
                                plot = "pdf", # 结果图片的导出类型，可以选择"png"或者"pdf" 
                                title = "Consensus Cluster")
Kvec <- 2:7
x1 = 0.1
x2 = 0.9
PAC <- rep(NA,length(Kvec))
names(PAC) <- paste("K =", Kvec, sep="")
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}
optK <- Kvec[which.min(PAC)]# 理想的K值
optK
#[1] 2


results[[3]]$consensusClass 

