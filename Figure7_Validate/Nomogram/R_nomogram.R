setwd("D:/MyProjects/scRNA_immune/Mac_Lung/Figure7_Validate/Nomogram")
rm(list = ls())

library(rms)
library(timeROC)
library(survival)
library(regplot)

risk<-read.table("D:/MyProjects/scRNA_immune/Mac_Lung/Figure5_M1_prognosis/risk.txt",
                 sep="\t",header=T, row.names = 1)
head(risk)
risk$sample_name<-row.names(risk)

for(i in 1:nrow(risk)){
  tmp <- strsplit(risk$sample_name[i], split = "_")[[1]][1]
  risk$sample_name[i] <- tmp
}

result_time<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/Figure6_prepareIndep/final_clinicaldata.csv",header=T)
head(result_time)
for(i in 1:nrow(result_time)){
  tmp <- strsplit(result_time$sample_name[i], split = "_")[[1]][1]
  result_time$sample_name[i] <- tmp
}

risk$sample_name<-row.names(risk)
risk$sample_name<-substr(risk$sample_name,1,12)
head(risk)


merged_data<-merge(risk, result_time, by = "sample_name")
merged_data<-subset(merged_data, select = c("sample_name","futime.x","fustat.x","ADAM19","ICAM3","WIPF1","LAP3","riskScore.x",
                                            "risk","gender","age","stage"))
colnames(merged_data)<-c("sample_name","futime","Status","ADAM19","ICAM3","WIPF1","LAP3","riskScore","risk","Gender",
                         "Age","Stage")

head(merged_data)
data_selected<-merged_data
colnames(data_selected)<-c("sample_name","futime","Status","ADAM19","ICAM3","WIPF1","LAP3","riskScore","risk","Gender",
                           "Age","Stage")
rt<-data_selected
head(rt)
for(i in 1:(dim(rt)[1])){
  rt$sample_name_new[i]<-paste(rt$sample_name[i],i,sep="_")
}
rownames(rt)<-rt$sample_name_new
rt<-subset(rt, select = c("futime", "Status", "Age", "Gender", "Stage", "riskScore"))
colnames(rt)<-c("futime", "fustat", "Age", "Gender", "Stage", "riskScore")

head(rt)
write.csv(rt, "nomogram_risk.csv", quote = F)

ddist <- datadist(rt)
options(datadist='ddist')

f_cph <- cph(Surv(futime,fustat) ~ Age+Gender+Stage+riskScore,
             x=T, y=T, surv=T,
             data=rt)
f_cph
#Frequencies of Missing Values Due to Each Variable
#Surv(futime, fustat)                  Age               Gender                Stage            riskScore 
#                   0                   16                    0                   15                    0 

#Cox Proportional Hazards Model

#cph(formula = Surv(futime, fustat) ~ Age + Gender + Stage + riskScore, 
    data = rt, x = T, y = T, surv = T)


#Model Tests     Discrimination    
#Indexes    
#Obs      1399    LR chi2    114.54     R2       0.079    
#Events    579    d.f.            4    R2(4,1399)0.076    
#Center 2.2898    Pr(> chi2) 0.0000     R2(4,579)0.174    
#Score chi2 134.93     Dxy      0.275    
#Pr(> chi2) 0.0000                       

#Coef   S.E.   Wald Z Pr(>|Z|)
#Age       0.0145 0.0046 3.13   0.0017  
#Gender    0.3198 0.0880 3.63   0.0003  
#Stage     0.3729 0.0450 8.28   <0.0001 
#riskScore 0.5100 0.0751 6.79   <0.0001 

#查看多因素Cox分析结果，最下方可见其对应的Coef、p值
print(f_cph)
ddist <- datadist(rt)
options(datadist='ddist')
med  <- Quantile(f_cph)
surv <- Survival(f_cph) 

pdf("nomogram_distribution.pdf")
regplot(f_cph,        #对观测2的六个指标在列线图上进行计分展示        
#        observation=rt[6,], #也可以不展示        
        points=TRUE,        
        plots=c("density","no plot"),        #预测1年和2年的死亡风险，此处单位是day        
        failtime = c(3,5,10),
        odds=F,        
        droplines=F,        
        leftlabel=T,        
        prfail = TRUE, #cox回归中需要TRUE        
        showP = T, #是否展示统计学差异        
        #droplines = F,#观测2示例计分是否画线        #    
#        colors = mycol, #用前面自己定义的颜色        
        rank="range", #根据统计学差异的显著性进行变量的排序        
        interval="confidence",        
        title="Cox regression") #展示观测的可信区间## [[1]]##   
dev.off()

pdf("nomogram.pdf",width=8, height=6)
plot(nomogram(f_cph, fun=list(function(x) surv(3, x),
                              function(x) surv(5, x),
                              function(x) surv(10, x)),
              funlabel=c("3-year Survival Probability", 
                         "5-year Survival Probability",
                         "10-year Survival Probability"))
)
dev.off()

pdf("Age_distribution.pdf")
theme_set(theme_bw())
p<-ggplot(rt,aes(x=Age,fill=T,alpha = 0.25))+geom_density()
p_bottom=p+theme(legend.position = "top")
p_bottom
dev.off()

pdf("Gender_distribution.pdf")
barplot(table(rt$Gender))
dev.off()

pdf("Stage_distribution.pdf")
barplot(table(rt$Stage))
dev.off()

pdf("RiskScore_distribution.pdf")
p<-ggplot(rt,aes(x=riskScore,fill=T,alpha = 0.25))+geom_density()
p_bottom=p+theme(legend.position = "top")
p_bottom
dev.off()

##############################################################################################
risk<-read.table("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/04_multivariate_Cox_regression_analysis/risk.txt",
                 sep="\t",header=T, row.names = 1)
head(risk)
result_time<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/Metastatic_nonMetastatic/result_time.csv",header=T)
head(result_time)

risk$sample_name<-row.names(risk)
risk$sample_name<-substr(risk$sample_name,1,12)
head(risk)
merged_data<-merge(risk, result_time, by.x = "sample_name", by.y = "submitter_id")
merged_data<-subset(merged_data, select = c("sample_name","futime.x","fustat.x","TRIP10","NGFR","SLC48A1","SRMS","riskScore",
                                            "risk","gender","age_at_index","ajcc_pathologic_stage"))
colnames(merged_data)<-c("sample_name","futime","Status","TRIP10","NGFR","SLC48A1","SRMS","riskScore","risk","Gender",
                         "Age","Stage")
merged_data$Gender_num<-ifelse(merged_data$Gender=="female","2",
                               ifelse(merged_data$Gender=="male","1",
                                      "NA"))
merged_data$Age_num<-ifelse(merged_data$Age>=60,"1",
                            ifelse(merged_data$Age<60,"2",
                                   "NA"))
merged_data$Stage_num<-ifelse(merged_data$Stage=="Stage IA", "1",
                              ifelse(merged_data$Stage=="Stage IIA", "2",
                                     ifelse(merged_data$Stage=="Stage IIB", "2",
                                            ifelse(merged_data$Stage=="Stage IIC", "2",
                                                   ifelse(merged_data$Stage=="Stage IIIA", "3",
                                                          ifelse(merged_data$Stage=="Stage IIIB", "3",
                                                                 ifelse(merged_data$Stage=="Stage IIIC", "3",
                                                                        ifelse(merged_data$Stage=="Stage IVA", "4",
                                                                               ifelse(merged_data$Stage=="Stage IVB", "4",
                                                                                      ifelse(merged_data$Stage=="Stage I", "1",
                                                                                             ifelse(merged_data$Stage=="Stage II", "2",
                                                                                                    ifelse(merged_data$Stage=="Stage III", "3",
                                                                                                           ifelse(merged_data$Stage=="Stage IV", "4",
                                                                                                                  "NA")))))))))))))

for(i in 1:(dim(merged_data)[1])){
  merged_data$sample_name_new[i]<-paste(merged_data$sample_name[i],i,sep="_")
}
rownames(merged_data)<-merged_data$sample_name_new
head(merged_data)
merged_data_selected<-subset(merged_data,select = c("futime", "Status", "Age_num", "Gender_num", "Stage_num", "riskScore"))
merged_data_selected$futime<-as.numeric(merged_data_selected$futime)
merged_data_selected$Status<-as.numeric(merged_data_selected$Status)
merged_data_selected$Age_num<-as.numeric(merged_data_selected$Age_num)
merged_data_selected$Gender_num<-as.numeric(merged_data_selected$Gender_num)
merged_data_selected$Stage_num<-as.numeric(merged_data_selected$Stage_num)
merged_data_selected_rmNA<-na.omit(merged_data_selected)
f<-coxph(Surv(futime,Status)~Age_num+Gender_num+Stage_num+riskScore,merged_data_selected_rmNA)
merged_data_selected_rmNA$pr_failure3<-c(1-(summary(survfit(f,newdata=merged_data_selected_rmNA),times=3)$surv))
merged_data_selected_rmNA$pr_failure5<-c(1-(summary(survfit(f,newdata=merged_data_selected_rmNA),times=5)$surv))
merged_data_selected_rmNA$pr_failure10<-c(1-(summary(survfit(f,newdata=merged_data_selected_rmNA),times=10)$surv))
head(merged_data_selected_rmNA)
library(dplyr)
library(dcurves)
pdf("DCA_year3.pdf")
dca(Surv(futime,Status)~pr_failure3,
    data = merged_data_selected_rmNA,
    time = 3,
    thresholds = 1:50/100) %>%
  plot(smooth = T)
dev.off()

pdf("DCA_year5.pdf")
dca(Surv(futime,Status)~pr_failure5,
    data = merged_data_selected_rmNA,
    time = 5,
    thresholds = 1:50/100) %>%
  plot(smooth = T)
dev.off()

pdf("DCA_year10.pdf")
dca(Surv(futime,Status)~pr_failure10,
    data = merged_data_selected_rmNA,
    time = 10,
    thresholds = 1:50/100) %>%
  plot(smooth = T)
dev.off()


