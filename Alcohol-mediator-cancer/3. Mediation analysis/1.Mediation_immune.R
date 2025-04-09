.libPaths("D:/PackagesR/R_LIBS")
rm(list=ls())


library(TwoSampleMR)
library(kableExtra)
library(ggplot2)
library(cowplot)
library(dplyr) 
library(tidyr)

load(file = "3. Mediation analysis/res_alcohol_cancer_final.Rdata")
immid <- paste0('ebi-a-GCST9000',c(1391:2121))
df<-res_ivw_pancancer[res_ivw_pancancer$pval<0.05,]

load(file = "res _cancertoalcohol.Rdata")
res_cancer_to_alcohol_ivw<-res_cancer_to_alcohol[res_cancer_to_alcohol$method=="Inverse variance weighted",]
res_median_total<-data.frame()
for (i in 1:nrow(df)) {
  row<-which(res_cancer_to_alcohol_ivw$id.exposure==df$id.outcome[i] & res_cancer_to_alcohol_ivw$id.outcome==df$id.exposure[i])
  res<-data.frame(
    X=df$id.exposure[i],
    Y=df$id.outcome[i],
    M=paste0('ebi-a-GCST9000',c(1391:2121)),
    beta_total=df$b[i],
    p_total=df$pval[i],
    p_reverse=res_cancer_to_alcohol_ivw$pval[row]
    )
  res_median_total<-rbind(res_median_total,res)
}

save(res_median_total,file = "alcohol to cancer total.Rdata")
res_median_total_p<-res_median_total[res_median_total$p_reverse>0.05,]

res_alcohol_cancer_median<-res_median_total_p%>%
  mutate( beta1=NA,
          p1=NA,
          beta2=NA,
          p2=NA)

test<-res_alcohol_cancer_median[res_alcohol_cancer_median$M=="ebi-a-GCST90001391",c(1,2)]
for (i in 1:nrow(test)) {
  a<-test$X[i]
  b<-test$Y[i]
  #X
  expo_rt3 <- extract_instruments(outcome = a, p1 = 5e-6,
                                  clump = T, p2 = 5e-8, 
                                  r2 = 0.001, kb = 10000)
  
  R<- TwoSampleMR::get_r_from_bsen(expo_rt3$beta.exposure,expo_rt3$se.exposure,expo_rt3$samplesize.exposure)
  if (!is.na(R)[1]) {
    F_value <- as.numeric(R*R)*(as.numeric(expo_rt3$samplesize.exposure) - 2)/(1 - as.numeric(R*R))
    expo_rt3<- cbind(expo_rt3, R2 = as.numeric(R*R), F_value = F_value)
    expo_rt3 <- expo_rt3[as.numeric(expo_rt3$F_value) > 10,]
  } 
  for (j in c(1:731)) {
    c<-res_alcohol_cancer_median$M[j]
    outc_rt3 <- tryCatch({  
      extract_outcome_data(expo_rt3$SNP, outcomes = c)  
    }, error = function(e) {  
      NA })  
    harm_rt3 <- tryCatch({  
      harmonise_data(
        exposure_dat =  expo_rt3, 
        outcome_dat = outc_rt3,action=2)
    }, error = function(e) {  
      NA })  
    mr_result3<- tryCatch({  
      mr(harm_rt3)    }, error = function(e) {  
      NA })  
    mr_result3_ivw<-tryCatch({  
      mr_result3[mr_result3$method=="Inverse variance weighted",]
    }, error = function(e) {  
      NA })  
    tryCatch({  
      res_alcohol_cancer_median$beta1[731*(i-1)+j]<-mr_result3_ivw$b
      res_alcohol_cancer_median$p1[731*(i-1)+j]<-mr_result3_ivw$pval
    }, error = function(e) {  
      NA })  
    ###### step 4 X' to Y TSMR(得到beta2) ######
    #M
    expo_rt4 <- tryCatch({  
      extract_instruments(outcome = c, p1 = 1e-5,
                          clump = T, p2 = 5e-8, 
                          r2 = 0.001, kb = 10000)
    }, error = function(e) {  
      NA })  
    R<- TwoSampleMR::get_r_from_bsen(expo_rt4$beta.exposure,expo_rt4$se.exposure,expo_rt4$samplesize.exposure)
    F_value <- as.numeric(R*R)*(as.numeric(expo_rt4$samplesize.exposure) - 2)/(1 - as.numeric(R*R))
    expo_rt4<- cbind(expo_rt4, R2 = as.numeric(R*R), F_value = F_value)
    expo_rt4 <- expo_rt4[as.numeric(expo_rt4$F_value) > 10,]
    
    outc_rt4 <- tryCatch({  
      extract_outcome_data(expo_rt4$SNP, outcomes =b)
    }, error = function(e) {  
      NA })  
    harm_rt4 <- tryCatch({  
      harmonise_data(
        exposure_dat =  expo_rt4, 
        outcome_dat = outc_rt4,action=2)
    }, error = function(e) {  
      NA })  
    mr_result4<- tryCatch({  
      mr(harm_rt4)
    }, error = function(e) {  
      NA })  
    mr_result4_ivw<-tryCatch({  
      mr_result4[mr_result4$method=="Inverse variance weighted",]
    }, error = function(e) {  
      NA })  
    tryCatch({  
      res_alcohol_cancer_median$beta2[731*(i-1)+j]<-mr_result4_ivw$b
      res_alcohol_cancer_median$p2[731*(i-1)+j]<-mr_result4_ivw$pval
    }, error = function(e) {  
      NA })  
   
    
  }
}
#addition<-res_alcohol_cancer_median
save(res_alcohol_cancer_median,file = "Mediation alcohol to cancer.Rdata")
res_alcohol_cancer_median_p<-res_alcohol_cancer_median %>%
  filter(p1 < 0.05 & p2 < 0.05)


median_res<-res_alcohol_cancer_median_p%>%
  mutate(se1 = NA) %>%
  mutate(se2 = NA) 

for (t in 1:nrow(median_res)) {
  a<-median_res$X[t]
  b<-median_res$Y[t]
  #X
  expo_rt3 <- extract_instruments(outcome = a, p1 = 5e-6,
                                  clump = T, p2 = 5e-8, 
                                  r2 = 0.001, kb = 10000)
  
  R<- TwoSampleMR::get_r_from_bsen(expo_rt3$beta.exposure,expo_rt3$se.exposure,expo_rt3$samplesize.exposure)
  if (!is.na(R)[1]) {
    F_value <- as.numeric(R*R)*(as.numeric(expo_rt3$samplesize.exposure) - 2)/(1 - as.numeric(R*R))
    expo_rt3<- cbind(expo_rt3, R2 = as.numeric(R*R), F_value = F_value)
    expo_rt3 <- expo_rt3[as.numeric(expo_rt3$F_value) > 10,]
  } 
  c<-median_res$M[t]
  outc_rt3 <- tryCatch({  
    extract_outcome_data(expo_rt3$SNP, outcomes = c)  
  }, error = function(e) {  
    NA })  
  harm_rt3 <- tryCatch({  
    harmonise_data(
      exposure_dat =  expo_rt3, 
      outcome_dat = outc_rt3,action=2)
  }, error = function(e) {  
    NA })  
  mr_result3<- tryCatch({  
    mr(harm_rt3)    }, error = function(e) {  
      NA })  
  mr_result3_ivw<-tryCatch({  
    mr_result3[mr_result3$method=="Inverse variance weighted",]
  }, error = function(e) {  
    NA })  
 median_res$se1<-mr_result3_ivw$se
  expo_rt4 <- tryCatch({  
    extract_instruments(outcome = c, p1 = 1e-5,
                        clump = T, p2 = 5e-8, 
                        r2 = 0.001, kb = 10000)
  }, error = function(e) {  
    NA })  
  R<- TwoSampleMR::get_r_from_bsen(expo_rt4$beta.exposure,expo_rt4$se.exposure,expo_rt4$samplesize.exposure)
  F_value <- as.numeric(R*R)*(as.numeric(expo_rt4$samplesize.exposure) - 2)/(1 - as.numeric(R*R))
  expo_rt4<- cbind(expo_rt4, R2 = as.numeric(R*R), F_value = F_value)
  expo_rt4 <- expo_rt4[as.numeric(expo_rt4$F_value) > 10,]
  
  outc_rt4 <- tryCatch({  
    extract_outcome_data(expo_rt4$SNP, outcomes =b)
  }, error = function(e) {  
    NA })  
  harm_rt4 <- tryCatch({  
    harmonise_data(
      exposure_dat =  expo_rt4, 
      outcome_dat = outc_rt4,action=2)
  }, error = function(e) {  
    NA })  
  mr_result4<- tryCatch({  
    mr(harm_rt4)
  }, error = function(e) {  
    NA })  
  mr_result4_ivw<-tryCatch({  
    mr_result4[mr_result4$method=="Inverse variance weighted",]
  }, error = function(e) {  
    NA })  
  median_res$se2[t]<-mr_result4_ivw$se
}

# two step 
###### step 1 X to Y TSMR
###### step 2 Y to X TSMR
###### step 3 X to X' TSMR
###### step 4 X' to Y TSMR

#1)Delta:
product_method_Delta <- function(EM_beta, EM_se, MO_beta, MO_se, verbose=F){
  # calculate indirect effect of exposure on outcome (via mediator) 
  # i.e. how much mediator accounts for total effect of exposure on outcome effect
  
  # this function can run either of two method, depending on what MO data df you supply 
  # calculate indirect effect beta
  EO <- EM_beta * MO_beta
  if (verbose) {print(paste("Indirect effect = ", round(EM_beta, 2)," x ", round(MO_beta,2), " = ", round(EO, 3)))}

  # Calculate CIs using RMediation package
  CIs = medci(EM_beta, MO_beta, EM_se, MO_se, type="dop")
  
  # put data into a tidy df
  df <-data.frame(b = EO,
                  se = CIs$SE,
                  lo_ci = CIs$`95% CI`[1],
                  up_ci= CIs$`95% CI`[2])
  # calculate OR
  df$or        <-  exp(df$b)
  df$or_lci95 <- exp(df$lo_ci)
  df$or_uci95 <- exp(df$up_ci)
  
  df<-round(df,3)
  return(df)
}
# method 1
# INDIRECT = TOTAL (exposure -> mediator) x TOTAL (mediator -> outcome)
product_method_Delta(EM_beta, EM_se, Mo_beta_total, Mo_se_total)
# method 2
# INDIRECT = TOTAL (exposure -> mediator) x DIRECT (of mediator , mvmr)
#product_method_Delta(EM_beta, EM_se, Mo_beta, Mo_se)

#2)误差传染法:
product_method_PoE <- function(EM_beta, EM_se, MO_beta, MO_se, verbose=F){
  # calculate indirect effect of exposure on outcome (via mediator) 
  # i.e. how much mediator accounts for total effect of exposure on outcome effect
  
  # this function can run either of two method, depending on what MO data df you supply 

  # calculate indirect effect beta
  EO_beta <- EM_beta * MO_beta
  
  if (verbose) {print(paste("Indirect effect = ", round(EM_beta, 2)," x ", round(MO_beta,2), " = ", round(EO, 3)))}
  
  
  # calculate SE of indirect effect 
  ### using propagation of errors method
  # SE of INDIRECT effect (difference) = sqrt(SE EM^2 + SE MO^2) 
  EO_se = round(sqrt(EM_se^2 + MO_se^2), 4)
  if (verbose) {print(paste("SE of indirect effect = sqrt(",
                            round(EM_se, 2),"^2 + ", round(MO_se,2), 
                            "^2) = ", indirect_se))}
  
  
  # put data into a tidy df
  df <-data.frame(b= EO_beta,
                  se = EO_se)
  
  # calculate CIs and OR
  df$lo_ci    <- df$b - 1.96 * df$se
  df$up_ci    <- df$b + 1.96 * df$se
  df$or        <-  exp(df$b)
  df$or_lci95 <- exp(df$lo_ci)
  df$or_uci95 <- exp(df$up_ci)
  
  df<-round(df,3)
  return(df)
}


median_alcohol_cancer<-median_res %>%
  mutate(indirect_effect=NA)%>%
  mutate(lo_ci=NA)%>%
  mutate(up_ci=NA)%>%
  mutate(or=NA)%>%
  mutate(or_lci95=NA)%>%
  mutate(or_uci95=NA)%>%
  mutate(direct_effect=NA)%>%
  mutate(mediated_proportion=NA)%>%
  mutate(p_value=NA)



for (t in 1:nrow(median_alcohol_cancer)) {
  df<-product_method_PoE(median_alcohol_cancer$beta1[t], median_alcohol_cancer$se1[t], 
                           median_alcohol_cancer$beta2[t], median_alcohol_cancer$se2[t])
  median_alcohol_cancer$indirect_effect[t]<-df$b
  median_alcohol_cancer$lo_ci[t]<-df$lo_ci
  median_alcohol_cancer$up_ci[t]<-df$up_ci
  median_alcohol_cancer$or[t]<-df$or
  median_alcohol_cancer$or_lci95[t]<-df$or_lci95
  median_alcohol_cancer$or_uci95[t]<-df$or_uci95
  median_alcohol_cancer$direct_effect[t]<-median_alcohol_cancer$beta_total[t]-df$b
  median_alcohol_cancer$mediated_proportion[t]<-df$b/median_alcohol_cancer$beta_total[t]
  Z<-df$b/df$se
  P<-2*pnorm(q=abs(Z), lower.tail=FALSE)
  median_alcohol_cancer$p_value[t]<-P
  
}

load(file = "ID_all.Rdata")
load(file = "gwasid_alcohol.Rdata")
load(file = "gwasid_cancer.Rdata")

data_immune<-idTotrait[idTotrait$id %in% immid,]

data1<-data_immune[,c(1:2)]
colnames(data1)<-c('M',"Immune_cell")
res1<-merge(data1,median_alcohol_cancer,by="M")
data2<-data_cancer[,c(1:2)]
colnames(data2)<-c('Y',"Cancer")
res2<-merge(data2,res1,by="Y")
data3<-data_alcohol[,c(1:2)]
colnames(data3)<-c('X',"Alcohol")
res_median_immune<-merge(data3,res2,by="X")

save(res_median_immune,file = "alcohol to immune to cancer.Rdata")
write.csv(res_median_immune,file = "alcohol to immune to cancer.csv")







#product_method_PoE(EM_beta, EM_se, Mo_beta_total, Mo_se_total)
df<-product_method_PoE(median_alcohol_cancer$beta1[t], median_alcohol_cancer$se1[t], 
                       median_alcohol_cancer$beta2[t], median_alcohol_cancer$se2[t])
# beta12=beta1*beta2
Indirect_effect <- result_or3$b*result_or4$b
#  beta_dir=beta_all-beta12
Direct_effect<-result_or$b - Indirect_effecteffect
#Z
Z=Indirect_effect/se
#p
P=2*pnorm(q=abs(Z), lower.tail=FALSE)

Mediated_Proportion<- (Indirect_effecteffect/result_or$b)
lci_p=lci/beta_all
uci_p=uci/beta_all
