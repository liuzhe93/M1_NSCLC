setwd("D:/MyProjects/scRNA_immune/Mac_Lung/Figure5_M1_prognosis")
rm(list=ls())

library("glmnet")
library("survival")

genes<-read.table("intersectGenes.txt")
exp<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ExprData/merged_expr_dat_remove_batch_effect.csv",
              header = T,row.names = 1)
exp_selected<-exp[genes$V1,]
exp_selected_t<-t(exp_selected)
exp_selected_t$sample_name<-rownames(exp_selected_t)
exp_selected_t<-as.data.frame(exp_selected_t)
exp_selected_t$sample_name<-rownames(exp_selected_t)
exp_selected_t$sample_name<-substr(exp_selected_t$sample_name,1,12)
exp_selected_t$sample_name<-gsub("\\.","-",exp_selected_t$sample_name)

result_time<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ClinicalData/merged_clinicalinfor.csv")
dim(result_time)
#[1] 1432    3
merged_re<-merge(exp_selected_t, result_time, by = "sample_name")
dim(merged_re)
#[1] 1490   14
merged_re_sort<-merged_re[order(merged_re[,13],decreasing=TRUE),]
merged_re_sort_rmNA<-na.omit(merged_re_sort)
dim(merged_re_sort_rmNA)
#[1] 1475   14
write.csv(merged_re_sort_rmNA, "GeneExp_withSurvival.csv", quote = F)

for(i in 1:nrow(merged_re_sort_rmNA)){
  merged_re_sort_rmNA$sample_name[i]<-paste(merged_re_sort_rmNA$sample_name[i],i,sep = "_")
}

row.names(merged_re_sort_rmNA)<-merged_re_sort_rmNA$sample_name
merged_re_sort_rmNA<-merged_re_sort_rmNA[,-1]

rt=merged_re_sort_rmNA

rt$futime[rt$futime<=0]=1

rt<-subset(rt, select = c("futime","fustat",colnames(rt)[1:(ncol(rt)-2)]))


x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))

fit <- glmnet(x, y, family = "cox", maxit = 1000)
pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene=c("futime","fustat",lassoGene)
lassoSigExp=rt[,lassoGene]
lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSigExp,file="lassoSigExp.txt",sep="\t",row.names=F,quote=F)
