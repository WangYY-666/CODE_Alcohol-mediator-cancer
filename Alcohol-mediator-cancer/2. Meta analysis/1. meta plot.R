
rm(list=ls())

# Loading necessary libraries for data manipulation and analysis
library(readr)      # For reading data files
library(vroom)      # For fast data reading
library(tidyr)      # For data tidying
library(tibble)     # For modern data frames
library(dplyr)      # For data manipulation
library(TwoSampleMR) # For Mendelian Randomization (MR) analysis
library(ggplot2)    # For data visualization
library(gt)         # For creating tables
library(meta)       # For meta-analysis
library(ComplexHeatmap) # For creating complex heatmaps
library(reshape2)   # For reshaping data
library(grid)       # For grid graphics

load(file = "2. Meta analysis/res _alcoholtocancer.Rdata")

res_ivw<-res_alcohol_to_cancer[res_alcohol_to_cancer$method=="Inverse variance weighted",]
res_group<-read.table(file="2. Meta analysis/alcohol group.txt",sep="\t",header = TRUE)
colnames(res_group)[2]<-"id.exposure"

#Taking BLCA as an example
res_ivw_BLCA<-res_ivw[grepl("Bladder cancer",res_ivw$outcome),]
test<-res_ivw_BLCA


test<-merge(res_group,test,by="id.exposure")
test<-test[test$group %in% c("Alcohol","Former alcohol drinker","Alcohol consumption","Alcohol abuse"),]
test$group<-factor(test$group, levels = c("Alcohol","Former alcohol drinker","Alcohol consumption","Alcohol abuse"))
test<-test[order(test$group),]
f<-function(x) {
  if(x<0.01)
    sprintf("%.2e", x)
  else
    round(x,2)
}
test[,c(8,9)]<-round(test[,c(8,9)],4)
test$pval<-sapply(test$pval, f)
test$combined <- paste(test$b, "(", test$pval, ")", sep = " ")
#有异质性(p<0.1 或 I^2>50%)应需用随机效应模型，否则可选用固定效应模型
meta1<-metagen(
  TE = log(test$or),
  seTE = (log(test$or_uci95)-log(test$or_lci95))/(2*1.96),
  #studlab = paste(test$id.exposure),  # Use only the method as the label
  studlab = paste(test$id.exposure,test$id.outcome,sep = "→"),
  #random = TRUE,common = FALSE,
  random = FALSE,common = TRUE,
  data = test,
  sm = "OR", #Specify the effect measure as Odds Ratio
  subgroup = group, print.subgroup.name = FALSE,
  test.subgroup = FALSE)
forest(meta1)

