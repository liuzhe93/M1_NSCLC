setwd("D:/MyProjects/scRNA_immune/Mac_Lung/Figure2_Mac_WGCNA/ClinicalData")
remove(list=ls())

###################################GSE13213 clinical data pre-operation#######################################
GSE13213 <- read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/external/GSE13213/clinical_GSE13213.csv",header = F, row.names = 1)
GSE13213_t <- t(GSE13213)
GSE13213_t <- as.data.frame(GSE13213_t)
GSE13213_t$Age <- gsub("Age: ","",GSE13213_t$Age)
GSE13213_t$Gender <- gsub("Sex: ","",GSE13213_t$Gender)
GSE13213_t$Histology <- gsub("Histology: ","",GSE13213_t$Histology)
GSE13213_t$Smoking_BI <- gsub("Smoking \\(BI\\): ","",GSE13213_t$Smoking_BI)
GSE13213_t$TNM_Pathological <-gsub("TNM \\(Pathological\\): ","",GSE13213_t$TNM_Pathological)
GSE13213_t$Stage_Pathological <- gsub("Stage \\(Pathological \\): ","",GSE13213_t$Stage_Pathological)
GSE13213_t$Status <- gsub("Status: ","",GSE13213_t$Status)
GSE13213_t$Survival_days <- gsub("Survival \\(days\\): ","",GSE13213_t$Survival_days)
GSE13213_t$Evidence_of_relapse <- gsub("Evidence of relapse: ","",GSE13213_t$Evidence_of_relapse)
GSE13213_t$EGFR_status <- gsub("EGFR status: ","",GSE13213_t$EGFR_status)
GSE13213_t$K_ras_Status<-gsub("K-ras Status: ","",GSE13213_t$K_ras_Status)
GSE13213_t$p53_Status<-gsub("p53 Status: ","",GSE13213_t$p53_Status)

#0=alive, 1=dead.
GSE13213_t$fustat <- ifelse(GSE13213_t$Status=="Alive",0,1)
GSE13213_t$futime <- as.numeric(GSE13213_t$Survival_days)/365
GSE13213_t<-subset(GSE13213_t, select=c("GSM", "fustat", "futime", "Source", "Age", "Gender", "Histology", "Smoking_BI", "TNM_Pathological",
                                         "Stage_Pathological", "Evidence_of_relapse", "EGFR_status","K_ras_Status", "p53_Status"))
write.csv(GSE13213_t,"result_time_GSE13213.csv",quote = F, row.names=F)

###################################GSE31210 clinical data pre-operation#######################################
GSE31210 <- read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/external/GSE31210/clinical_GSE31210.csv",header = F, row.names = 1)
GSE31210_t <- t(GSE31210)
GSE31210_t <- as.data.frame(GSE31210_t)
GSE31210_t$Tissue <- gsub("tissue: ","",GSE31210_t$Tissue)
GSE31210_t$Age <- gsub("age \\(years\\): ","",GSE31210_t$Age)
GSE31210_t$Age <- gsub("age: ","",GSE31210_t$Age)
GSE31210_t$Gender <- gsub("gender: ","",GSE31210_t$Gender)
GSE31210_t$Smoking_status <- gsub("smoking status: ","",GSE31210_t$Smoking_status)
GSE31210_t$BI <- gsub("bi: ","",GSE31210_t$BI)
GSE31210_t$Pathological_stage <- gsub("pathological stage: ","",GSE31210_t$Pathological_stage)
GSE31210_t$Pstage_iorii <- gsub("pstage iorii: ","",GSE31210_t$Pstage_iorii)
GSE31210_t$Gene_alteration_status <- gsub("gene alteration status: ","",GSE31210_t$Gene_alteration_status)
GSE31210_t$Myc <- gsub("myc: ","",GSE31210_t$Myc)
GSE31210_t$Myc_copy <- gsub("myc_copy: ","",GSE31210_t$Myc_copy)
GSE31210_t$Cluster <- gsub("cluster: ","",GSE31210_t$Cluster)
GSE31210_t$Death <- gsub("death: ","",GSE31210_t$Death)
GSE31210_t$Time <- gsub("days before death\\/censor: ","",GSE31210_t$Time)
GSE31210_t<-GSE31210_t[1:226,]
#0=alive, 1=dead.
GSE31210_t$fustat <- ifelse(GSE31210_t$Death=="alive",0,1)
GSE31210_t$futime <- as.numeric(GSE31210_t$Time)/365

GSE31210_t<-subset(GSE31210_t, select=c("GSM", "fustat", "futime", "Tissue", "Age", "Gender", "Smoking_status","BI",
                                        "Pathological_stage", "Pstage_iorii", "Gene_alteration_status", "Myc","Myc_copy", "Cluster"))
write.csv(GSE31210_t,"result_time_GSE31210.csv",quote = F, row.names=F)

###################################LUAD clinical data pre-operation#######################################
surv_data<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/TCGA/LUAD/TCGAbiolinks-LUAD-clinical.csv")
surv_data<-subset(surv_data,select = c("submitter_id","ajcc_pathologic_stage","ajcc_pathologic_t","ajcc_pathologic_n",
                                       "ajcc_pathologic_m","gender","vital_status","age_at_index","cigarettes_per_day",
                                       "pack_years_smoked","days_to_last_follow_up","days_to_death"))
#cigarettes_per_day: The average number of cigarettes smoked per day.
#pack_years_smoked: Numeric computed value to represent lifetime tobacco exposure defined as number of cigarettes smoked per day x number of years smoked divided by 20.

surv_data$os<-ifelse(surv_data$vital_status=='Alive',surv_data$days_to_last_follow_up,surv_data$days_to_death)
library(dplyr)
surv_selected<-subset(surv_data, select=c("submitter_id","os","vital_status","gender","age_at_index","ajcc_pathologic_stage",
                                          "ajcc_pathologic_t","ajcc_pathologic_n","ajcc_pathologic_m","cigarettes_per_day",
                                          "pack_years_smoked"))
#0=alive, 1=dead.
surv_selected$fustat<-ifelse(surv_selected$vital_status=='Alive',0,1)
surv_selected$futime<-surv_selected$os/365
result_time<-subset(surv_selected, select=c("submitter_id","futime","fustat","gender","age_at_index","ajcc_pathologic_stage",
                                            "ajcc_pathologic_t","ajcc_pathologic_n","ajcc_pathologic_m","cigarettes_per_day",
                                            "pack_years_smoked"))
write.csv(result_time,"result_time_LUAD.csv",row.names=F)

###################################LUSC clinical data pre-operation#######################################
surv_data<-read.csv("D:/MyProjects/scRNA_immune/Mac_Lung/DataCollection/TCGA/LUSC/TCGAbiolinks-LUSC-clinical.csv")
surv_data<-subset(surv_data,select = c("submitter_id","ajcc_pathologic_stage","ajcc_pathologic_t","ajcc_pathologic_n",
                                       "ajcc_pathologic_m","gender","vital_status","age_at_index","cigarettes_per_day",
                                       "pack_years_smoked","days_to_last_follow_up","days_to_death"))
#cigarettes_per_day: The average number of cigarettes smoked per day.
#pack_years_smoked: Numeric computed value to represent lifetime tobacco exposure defined as number of cigarettes smoked per day x number of years smoked divided by 20.

surv_data$os<-ifelse(surv_data$vital_status=='Alive',surv_data$days_to_last_follow_up,surv_data$days_to_death)
library(dplyr)
surv_selected<-subset(surv_data, select=c("submitter_id","os","vital_status","gender","age_at_index","ajcc_pathologic_stage",
                                          "ajcc_pathologic_t","ajcc_pathologic_n","ajcc_pathologic_m","cigarettes_per_day",
                                          "pack_years_smoked"))
#0=alive, 1=dead.
surv_selected$fustat<-ifelse(surv_selected$vital_status=='Alive',0,1)
surv_selected$futime<-surv_selected$os/365
result_time<-subset(surv_selected, select=c("submitter_id","futime","fustat","gender","age_at_index","ajcc_pathologic_stage",
                                            "ajcc_pathologic_t","ajcc_pathologic_n","ajcc_pathologic_m","cigarettes_per_day",
                                            "pack_years_smoked"))
write.csv(result_time,"result_time_LUSC.csv",row.names=F)


#########################merge all clinical information###################################
clinical_GSE13213<-read.csv("result_time_GSE13213.csv",header=T)
clinical_GSE31210<-read.csv("result_time_GSE31210.csv",header=T)
TCGA_LUAD<-read.csv("result_time_LUAD.csv",header=T)
TCGA_LUSC<-read.csv("result_time_LUSC.csv",header=T)
clinical_GSE13213<-clinical_GSE13213[,1:3]
colnames(clinical_GSE13213)<-c("sample_name","fustat","futime" )
clinical_GSE31210<-clinical_GSE31210[,1:3]
colnames(clinical_GSE31210)<-c("sample_name","fustat","futime" )
TCGA_LUAD<-TCGA_LUAD[,1:3]
colnames(TCGA_LUAD)<-c("sample_name","futime","fustat" )
TCGA_LUSC<-TCGA_LUSC[,1:3]
colnames(TCGA_LUSC)<-c("sample_name","futime","fustat" )

merged_clinical<-rbind(clinical_GSE13213,clinical_GSE31210,TCGA_LUAD,TCGA_LUSC)
dim(merged_clinical)
#[1] 1432    3
write.csv(merged_clinical, "merged_clinicalinfor.csv", quote = F, row.names = F)





