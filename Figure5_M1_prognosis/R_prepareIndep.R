setwd("D:/MyProjects/scRNA_immune/Mac_Lung/Figure5_M1_prognosis")
rm(list=ls())

clinical_data<-read.table("clinical_rmNA.txt",header=T)
risk_data<-read.table("risk.txt",header=T)
colnames(clinical_data)<-c("id", "futime", "fustat", "age", "gender", "stage", "T", "M", "N")
merge_data<-merge(clinical_data,risk_data,by="id")
merge_data<-merge_data %>%
  select(id,futime.x,fustat.x,age,gender,stage,T,M,N,riskScore)		
colnames(merge_data)<-c("id","futime","fustat","age","gender","stage","T","M","N","riskScore")
write.table(merge_data,"indepInput.txt",sep="\t",quote=F,row.names=F)


