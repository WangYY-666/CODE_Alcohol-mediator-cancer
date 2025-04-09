
rm(list=ls())

library(MRInstruments)
library(plyr)
library(dplyr)
library(TwoSampleMR)
library(writexl)

idTotrait <- available_outcomes()
IEUid<-as.data.frame(idTotrait)
rownames(IEUid)<-IEUid[,1]
load(file = "1. Bidirectional  MR/gwasid_alcohol.Rdata")
load(file = "1. Bidirectional  MR/gwasid_cancer.Rdata")

res_cancer_to_alcohol<-data.frame()
for (i in data_cancer$id) {
  expo_data <- extract_instruments(outcome = i, p1 = 5e-6,
                                   clump = T, p2 = 5e-8,
                                   r2 = 0.001, kb = 10000)
  #计算F检验值：方法1
  R<- TwoSampleMR::get_r_from_bsen(expo_data$beta.exposure,expo_data$se.exposure,expo_data$samplesize.exposure)
  F_value <- as.numeric(R*R)*(as.numeric(expo_data$samplesize.exposure) - 2)/(1 - as.numeric(R*R))
  expo_data<- cbind(expo_data, R2 = as.numeric(R*R), F_value = F_value)
  expo_data <- expo_data[as.numeric(expo_data$F_value) > 10,]
  for(j in data_alcohol$id[1:nrow(data_cancer)]){
    # 提取结果因变量数据
    outc_data <- extract_outcome_data(snps = expo_data$SNP, outcomes = j)
    # 调和暴露因素和结果数据
    harm_data <- harmonise_data(exposure_dat = expo_data,
                                outcome_dat = outc_data,
                                action = 2)
    # 进行MR分析
    mr_result <- mr(harm_data)
    # 生成优比值
    result_or <- generate_odds_ratios(mr_result)
    # 进行多重比较校正
    pleiotropy <- mr_pleiotropy_test(harm_data)
    #write.table(pleiotropy, file = paste0(filename,"/pleiotropy.txt"),sep = "\t", quote = F)
    # 进行异质性检验
    heterogeneity <- mr_heterogeneity(harm_data)
    heterQ<-data.frame("heterogeneity "=c(heterogeneity$Q_pval[1],"",heterogeneity$Q_pval[2],"",""))
    #write.table(heterogeneity, file = paste0(filename,"/heterogeneity.txt"),sep = "\t", quote = F)
    if(nrow(result_or)==1){
      res<-result_or%>%
        mutate(
          pleiotropy = NA,  
          heterogeneity. = NA                
        )
    }else{
      res<-cbind(result_or,"pleiotropy"=pleiotropy$pval,"heterogeneity"=heterQ)
    }
    res_cancer_to_alcohol<-rbind(res_cancer_to_alcohol,res)
  }
}


res_alcohol_to_cancer<-data.frame()
for (i in data_alcohol$id) {
  expo_data <- extract_instruments(outcome = i, p1 = 5e-6,
                                   clump = T, p2 = 5e-8,
                                   r2 = 0.001, kb = 10000)
  #计算F检验值：方法1
  R<- TwoSampleMR::get_r_from_bsen(expo_data$beta.exposure,expo_data$se.exposure,expo_data$samplesize.exposure)
  F_value <- as.numeric(R*R)*(as.numeric(expo_data$samplesize.exposure) - 2)/(1 - as.numeric(R*R))
  expo_data<- cbind(expo_data, R2 = as.numeric(R*R), F_value = F_value)
  expo_data <- expo_data[as.numeric(expo_data$F_value) > 10,]
  for(j in data_cancer$id[1:nrow(data_cancer)]){
    # 提取结果因变量数据
    outc_data <- extract_outcome_data(snps = expo_data$SNP, outcomes = j)
    # 调和暴露因素和结果数据
    harm_data <- harmonise_data(exposure_dat = expo_data,
                                outcome_dat = outc_data,
                                action = 2)
    # 进行MR分析
    mr_result <- mr(harm_data)
    # 生成优比值
    result_or <- generate_odds_ratios(mr_result)
    # 进行多重比较校正
    pleiotropy <- mr_pleiotropy_test(harm_data)
    #write.table(pleiotropy, file = paste0(filename,"/pleiotropy.txt"),sep = "\t", quote = F)
    # 进行异质性检验
    heterogeneity <- mr_heterogeneity(harm_data)
    heterQ<-data.frame("heterogeneity "=c(heterogeneity$Q_pval[1],"",heterogeneity$Q_pval[2],"",""))
    #write.table(heterogeneity, file = paste0(filename,"/heterogeneity.txt"),sep = "\t", quote = F)
    if(nrow(result_or)==1){
      res<-result_or%>%
        mutate(
          pleiotropy = NA,  
          heterogeneity. = NA                
        )
    }else{
      res<-cbind(result_or,"pleiotropy"=pleiotropy$pval,"heterogeneity"=heterQ)
    }
    res_alcohol_to_cancer<-rbind(res_alcohol_to_cancer,res)
  }
}

