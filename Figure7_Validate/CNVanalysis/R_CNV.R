setwd("D:/MyProjects/scRNA_immune/Mac_Lung/Figure7_Validate/CNVanalysis")
rm(list = ls())

library("ggplot2")
# https://www.cbioportal.org/study/summary?id=nsclc_tcga_broad_2016
# https://www.cbioportal.org/
# download CNA_Genes.txt file
cna_data<-read.table("CNA_Genes.txt",sep="\t",header = T)
geneList <- c("ADAM19", "ICAM3", "WIPF1", "LAP3")


gene_index<-t(t(geneList))
colnames(gene_index)<-"Gene"
gene_index<-as.data.frame(gene_index)
data_export<-merge(cna_data, gene_index, by = "Gene")
amp_freq<-subset(data_export,CNA=="AMP")
summary(as.numeric(sub("%","",amp_freq$Freq))/100)
del_freq<-subset(data_export,CNA=="HOMDEL")
summary(as.numeric(sub("%","",del_freq$Freq))/100)
del_freq$Freq<-paste0("-",del_freq$Freq)
amp_freq$Freq<-as.numeric(sub("%","",amp_freq$Freq))
del_freq$Freq<-as.numeric(sub("%","",del_freq$Freq))
data<-rbind(amp_freq,del_freq)

pdf("CNA.pdf",width=8,height = 6)
ggplot(data = data) +
  geom_col(aes(x = factor(Gene), y = Freq, fill = CNA)) +
  scale_y_continuous(breaks = seq(from = -8, to = 16,by = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = .5))
dev.off()

