setwd("D:/MyProjects/scRNA_immune/Mac_Lung/Figure5_M1_prognosis")
rm(list=ls())

library(survival)
library("survminer")
rt=read.table("risk.txt",header=T,sep="\t")

outFile="survival.pdf"

#获取最优cutoff
res.cut=surv_cutpoint(rt, time = "futime", event = "fustat",variables =c("riskScore"))
res.cut
res.cat=surv_categorize(res.cut)
fit=survfit(Surv(futime, fustat) ~riskScore, data = res.cat)

#比较高低表达生存差异
diff=survdiff(Surv(futime, fustat) ~riskScore,data =res.cat)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}

#绘制
surPlot=ggsurvplot(fit, 
                   data=res.cat,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title="RiskScore",
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf(file=outFile,onefile = FALSE,width = 8,height =5)
print(surPlot)
dev.off()


res.cut
#          cutpoint statistic
#riskScore 1.145914  5.236254
rt$group=ifelse(rt$riskScore>1.145914,"high","low")
write.table(rt, "risk_optimalCutoff.txt", sep="\t", quote=F, row.names=F)

