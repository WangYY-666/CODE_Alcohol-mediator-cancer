
rm(list=ls())

library(TwoSampleMR)
library(dplyr)
load(file = "gwasid_alcohol.Rdata")
load(file = "gwasid_cancer.Rdata")

res_steiger_test<-data.frame()
for (i in data_alcohol$id) {
  expo_data <- extract_instruments(outcome = i, p1 = 5e-6,
                                   clump = T, p2 = 5e-8,
                                   r2 = 0.001, kb = 10000)
  #计算F检验值
  if(is.na(expo_data$samplesize.exposure[1])) {
    expo_data$samplesize.exposure<- (data_alcohol[data_alcohol$id==i,]$ncase+data_alcohol[data_alcohol$id==i,]$ncontrol)}
  R<- TwoSampleMR::get_r_from_bsen(expo_data$beta.exposure,expo_data$se.exposure,expo_data$samplesize.exposure)
  F_value <- as.numeric(R*R)*(as.numeric(expo_data$samplesize.exposure) - 2)/(1 - as.numeric(R*R))
  expo_data<- cbind(expo_data, r.exposure = as.numeric(R), F_value = F_value)
  expo_data <- expo_data[as.numeric(expo_data$F_value) > 10,]
  for (j in data_cancer$id) {
    outc_data <- extract_outcome_data(snps = expo_data$SNP, outcomes = j)
    if(is.na(outc_data$samplesize.outcome[1])) {
      outc_data$samplesize.outcome<- (data_cancer[data_cancer$id==j,]$ncase+data_cancer[data_cancer$id==j,]$ncontrol)}
    if(outc_data$samplesize.outcome[1]== 0){
      outc_data$samplesize.outcome<- data_cancer[data_cancer$id==j,]$sample_size
    }
     harm_data <- harmonise_data(exposure_dat = expo_data,
                                outcome_dat = outc_data,
                                action = 2)
    
    harm_data$r.outcome = TwoSampleMR::get_r_from_bsen(outc_data$beta.outcome,outc_data$se.outcome,outc_data$samplesize.outcome)      # 结局的R
    steiger_test <- directionality_test(harm_data)
    res_steiger_test<-rbind(res_steiger_test,steiger_test)
  }
  
}


