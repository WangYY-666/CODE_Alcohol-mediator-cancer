
rm(list=ls())


library(ggplot2)
library("survival") #生存分析
library("survminer") #结果可视化
library(dplyr)

load(file = "Clinic data(ESCA). Rdata")
clin_ESCA<-clin_ESCA[clin_ESCA$group=="cancer",]
clin_ESCA$OS.time_month<-as.numeric(clin_ESCA$OS.time)/30
clin_ESCA<-clin_ESCA[substr(clin_ESCA$sample,14,16)=="01A",]
fit1 <- survfit(Surv(OS.time_month, OS) ~ alcohol_history, data = clin_ESCA)

load(file = "clinic data merged(ESCA).Rdata")
clin_ESCA<-clinc_ESCA
clin_ESCA<-clin_ESCA[clin_ESCA$group=="cancer",]
clin_ESCA$OS.time_month<-as.numeric(clin_ESCA$OS.time)/30
colnames(clin_ESCA)
#clin_ESCA<-clin_ESCA[clin_ESCA$history_of_esophageal_cancer=="NO",]
#clin_ESCA<-clin_ESCA[clin_ESCA$h_pylori_infection=="Never",]
clin_ESCA<-clin_ESCA[clin_ESCA$histological_type=="Esophagus Squamous Cell Carcinoma",]
clin_ESCA<-clin_ESCA[clin_ESCA$histological_type=="Esophagus Adenocarcinoma, NOS",]
#clin_ESCA<-clin_ESCA[clin_ESCA$treatment_prior_to_surgery=="No Treatment",]
fit1 <- survfit(Surv(OS.time_month, OS) ~ alcohol_history, data = clin_ESCA)
#------------------------------------------------------------------------------

#fit2<-survdiff(Surv(OS.time_month, OS) ~ alcohol_history, data = clin_ESCA)
ggsurvplot(fit1,
           pval = TRUE, 
           conf.int = TRUE,
           risk.table = TRUE, # Add risk table
          # risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(),
           xlab = "Time since diagnosis (month)", legend = c(.9, .8),
           legend.labs = c("No", "Yes"),legend.title = "Alcohol history",
           palette = c("skyblue", "red1"))

clin_CA<-clin_ESCA[clin_ESCA$histological_type=="Esophagus Adenocarcinoma, NOS",]
clin_SCC<-clin_ESCA[clin_ESCA$histological_type=="Esophagus Squamous Cell Carcinoma",]

save(clin_CA,clin_SCC,file = "clinic_data(types).Rdata")






