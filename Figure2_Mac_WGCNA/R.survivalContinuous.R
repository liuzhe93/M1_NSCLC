setwd("/Users/liuzhe/Desktop/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/SurvivalAnalysis/Continuous_median/")
remove(list=ls())

rt<-read.csv("/Users/liuzhe/Desktop/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv")
library("survival")
library("survminer")

############"B.cells.naive"#######################################
rt<-read.csv("/Users/liuzhe/Desktop/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv")
var="B.cells.naive"
rt=rt[,c("futime","fustat",var)]
group=ifelse(rt[,3]>median(rt[,3]),"High","Low")
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = rt) #剩余数目

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf("survival_B.cells.naive.pdf",onefile = FALSE,width=10)
print(surPlot)
dev.off()

############"B.cells.memory"######################################
rt<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv")
var="B.cells.memory"
rt=rt[,c("futime","fustat",var)]
group=ifelse(rt[,3]>median(rt[,3]),"High","Low")
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = rt) #剩余数目

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf("survival_B.cells.memory.pdf",onefile = FALSE,width=10)
print(surPlot)
dev.off()

############"Plasma.cells"########################################
rt<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv")
var="Plasma.cells"
rt=rt[,c("futime","fustat",var)]
group=ifelse(rt[,3]>median(rt[,3]),"High","Low")
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = rt) #剩余数目

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf("survival_Plasma.cells.pdf",onefile = FALSE,width=10)
print(surPlot)
dev.off()

############"T.cells.CD4.naive"###################################
rt<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv")
var="T.cells.CD4.naive"
rt=rt[,c("futime","fustat",var)]
group=ifelse(rt[,3]>median(rt[,3]),"High","Low")
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = rt) #剩余数目

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf("survival_T.cells.CD4.naive.pdf",onefile = FALSE,width=10)
print(surPlot)
dev.off()

############"T.cells.follicular.helper"###########################
rt<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv")
var="T.cells.follicular.helper"
rt=rt[,c("futime","fustat",var)]
group=ifelse(rt[,3]>median(rt[,3]),"High","Low")
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = rt) #剩余数目

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf("survival_T.cells.follicular.helper.pdf",onefile = FALSE,width=10)
print(surPlot)
dev.off()

############"T.cells.regulatory..Tregs."##########################
rt<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv")
var="T.cells.regulatory..Tregs."
rt=rt[,c("futime","fustat",var)]
group=ifelse(rt[,3]>median(rt[,3]),"High","Low")
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = rt) #剩余数目

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf("survival_T.cells.regulatory..Tregs..pdf",onefile = FALSE,width=10)
print(surPlot)
dev.off()

############"T.cells.gamma.delta"#################################
rt<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv")
var="T.cells.gamma.delta"
rt=rt[,c("futime","fustat",var)]
group=ifelse(rt[,3]>median(rt[,3]),"High","Low")
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = rt) #剩余数目

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf("survival_T.cells.gamma.delta.pdf",onefile = FALSE,width=10)
print(surPlot)
dev.off()

############"NK.cells.resting"####################################
rt<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv")
var="NK.cells.resting"
rt=rt[,c("futime","fustat",var)]
group=ifelse(rt[,3]>median(rt[,3]),"High","Low")
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = rt) #剩余数目

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf("survival_NK.cells.resting.pdf",onefile = FALSE,width=10)
print(surPlot)
dev.off()

############"Monocytes"###########################################
rt<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv")
var="Monocytes"
rt=rt[,c("futime","fustat",var)]
group=ifelse(rt[,3]>median(rt[,3]),"High","Low")
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = rt) #剩余数目

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf("survival_Monocytes.pdf",onefile = FALSE,width=10)
print(surPlot)
dev.off()

############"Macrophages.M0"######################################
rt<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv")
var="Macrophages.M0"
rt=rt[,c("futime","fustat",var)]
group=ifelse(rt[,3]>median(rt[,3]),"High","Low")
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = rt) #剩余数目

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf("survival_Macrophages.M0.pdf",onefile = FALSE,width=10)
print(surPlot)
dev.off()
############"Macrophages.M1"######################################
rt<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv")
var="Macrophages.M1"
rt=rt[,c("futime","fustat",var)]
group=ifelse(rt[,3]>median(rt[,3]),"High","Low")
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = rt) #剩余数目

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf("survival_Macrophages.M1.pdf",onefile = FALSE,width=10)
print(surPlot)
dev.off()
############"Macrophages.M2"######################################
rt<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv")
var="Macrophages.M2"
rt=rt[,c("futime","fustat",var)]
group=ifelse(rt[,3]>median(rt[,3]),"High","Low")
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = rt) #剩余数目

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf("survival_Macrophages.M2.pdf",onefile = FALSE,width=10)
print(surPlot)
dev.off()
############"Dendritic.cells.resting"#############################
rt<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv")
var="Dendritic.cells.resting"
rt=rt[,c("futime","fustat",var)]
group=ifelse(rt[,3]>median(rt[,3]),"High","Low")
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = rt) #剩余数目

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf("survival_Dendritic.cells.resting.pdf",onefile = FALSE,width=10)
print(surPlot)
dev.off()

############"Dendritic.cells.activated"###########################
rt<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv")
var="Dendritic.cells.activated"
rt=rt[,c("futime","fustat",var)]
group=ifelse(rt[,3]>median(rt[,3]),"High","Low")
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = rt) #剩余数目

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf("survival_Dendritic.cells.activated.pdf",onefile = FALSE,width=10)
print(surPlot)
dev.off()

############"Mast.cells.resting" #################################
rt<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv")
var="Mast.cells.resting"
rt=rt[,c("futime","fustat",var)]
group=ifelse(rt[,3]>median(rt[,3]),"High","Low")
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = rt) #剩余数目

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf("survival_Mast.cells.resting.pdf",onefile = FALSE,width=10)
print(surPlot)
dev.off()

############"Mast.cells.activated"################################
rt<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv")
var="Mast.cells.activated"
rt=rt[,c("futime","fustat",var)]
group=ifelse(rt[,3]>median(rt[,3]),"High","Low")
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = rt) #剩余数目

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf("survival_Mast.cells.activated.pdf",onefile = FALSE,width=10)
print(surPlot)
dev.off()

############"Eosinophils"#########################################
rt<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ImmuneInfiltration/ImmuneScore_withSurvival.csv")
var="Eosinophils"
rt=rt[,c("futime","fustat",var)]
group=ifelse(rt[,3]>median(rt[,3]),"High","Low")
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = rt) #剩余数目

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf("survival_Eosinophils.pdf",onefile = FALSE,width=10)
print(surPlot)
dev.off()
