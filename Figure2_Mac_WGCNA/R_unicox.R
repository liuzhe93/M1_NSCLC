setwd("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/UniCoxRegression")
remove(list=ls())

library(survival)
library(survminer)

#clrs <- fpColors(box="black",line="black", summary="black")             #定义森林图颜色

rt=read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv",header=T,row.names =1)
rt=subset(rt, select = c("sample_name", "futime", "fustat", "B.cells.naive", "B.cells.memory", "Plasma.cells", "T.cells.CD8",
                         "T.cells.CD4.naive","T.cells.CD4.memory.resting", "T.cells.CD4.memory.activated", "T.cells.follicular.helper",
                         "T.cells.regulatory..Tregs.", "T.cells.gamma.delta", "NK.cells.resting", "NK.cells.activated",
                         "Monocytes", "Macrophages.M0", "Macrophages.M1", "Macrophages.M2", "Dendritic.cells.resting", 
                         "Dendritic.cells.activated", "Mast.cells.resting","Mast.cells.activated", "Eosinophils", "Neutrophils"))
head(rt)
for(i in 1:(dim(rt)[1])){
  rt$sample_name_new[i]<-paste(rt$sample_name[i],i,sep="_")
}
rownames(rt)<-rt$sample_name_new
rt<-subset(rt, select = c("futime", "fustat", "B.cells.naive", "B.cells.memory", "Plasma.cells", "T.cells.CD8","T.cells.CD4.naive",
                          "T.cells.CD4.memory.resting", "T.cells.CD4.memory.activated", "T.cells.follicular.helper",
                          "T.cells.regulatory..Tregs.", "T.cells.gamma.delta", "NK.cells.resting", "NK.cells.activated",
                          "Monocytes", "Macrophages.M0", "Macrophages.M1", "Macrophages.M2", "Dendritic.cells.resting", 
                          "Dendritic.cells.activated", "Mast.cells.resting","Mast.cells.activated", "Eosinophils", "Neutrophils"))
head(rt)

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(outTab,file="uniCox.xls",sep="\t",row.names=F,quote=F)
