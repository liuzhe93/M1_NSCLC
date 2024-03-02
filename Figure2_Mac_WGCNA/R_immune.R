setwd("/Users/liuzhe/Desktop/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/")
remove(list=ls())

exp<-read.csv("../ExprData/merged_expr_dat_remove_batch_effect.csv",header = T)
write.table(exp, "exp.txt", sep = "\t", quote = F, row.names = F)

#################################immune cell infiltration analysis######################################
source("CIBERSORT.R")
if(T){
  TME.results = CIBERSORT("LM22.txt", 
                          "exp.txt" , 
                          perm = 1000, 
                          QN = F)
  save(TME.results,file = "ciber_NSCLC.Rdata")
}

#if(T){
#  TME.results = CIBERSORT("LM22.txt", 
#                          "D:/MyProjects/CTC/CTC_Prognosis_COAD/05_ImmuneLandscape/03_ImmuneInfiltration/exp.txt" , 
#                          perm = 1000, 
#                          QN = F)
#  save(TME.results,file = "ciber_NSCLC.Rdata")
#}


load("ciber_NSCLC.Rdata")
TME.results[1:4,1:4]
dim(TME.results)
#[1] 1510   25
re <- TME.results[,-(23:25)]
library(pheatmap)
re_df<-as.data.frame(re)
re_df$sample_name<-row.names(re_df)
re_df$sample_name<-substr(re_df$sample_name,1,12)
re_df$sample_name<-gsub("\\.","-",re_df$sample_name)
dim(re_df)
#[1] 1510   23

#######################for prognostic analysis, we should remove normal samples#########################

re_df_rm_nor<-rbind(re_df[1:117,],re_df[118:343,],re_df[364:900,],re_df[960:1461,])

#################################read clinical information data#########################################
result_time<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ClinicalData/merged_clinicalinfor.csv")
dim(result_time)
#[1] 1432    3
merged_re<-merge(re_df_rm_nor, result_time, by = "sample_name")
dim(merged_re)
#[1] 1490   25
merged_re_sort<-merged_re[order(merged_re[,25],decreasing=TRUE),]
merged_re_sort_rmNA<-na.omit(merged_re_sort)
dim(merged_re_sort_rmNA)
#[1] 1382   25
write.csv(merged_re_sort_rmNA, "ImmuneScore_withSurvival.csv", quote = F)




