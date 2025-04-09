
rm(list=ls())

load(file = "ESCA_tpm_clin (176 samples).Rdata")
load(file = "clinic_data(types).Rdata")
rt1<-as.data.frame(rt_tpm[,c(1,14:ncol(rt_tpm))])
rt<-rt1[rt1$sample %in% clin_CA[,2],]
rownames(rt)<-rt$sample
rt<-rt[,-1]


library("survival")
library("survminer")

library(IOBR)
library(tidyverse)
library(clusterProfiler)

library(dplyr)
exp<-rt_tpm[,-c(1:13)]
rownames(exp)<-rt_tpm[,1]
exp<-as.data.frame(t(exp))

#exp<-exp[,colnames(exp) %in% clin_SCC$sample]
dim(exp)
tme_deconvolution_methods
#exp<-log2(exp+1)

# MCPcounter
mcpcounter <- deconvo_tme(eset = exp,tumor = TRUE,
                             method = "mcpcounter"
)
epic <- deconvo_tme(eset = exp,
                       method = "epic",
                       tumor = TRUE,
                       arrays = F
)
# xcell
xcell <- deconvo_tme(eset = exp,
                        method = "xcell",
                        arrays = F
)
# CIBERSORT
#perm：表示置换次数，数字越大运行时间越长，一般文章都设置为1000；
cibersort_tpm <- deconvo_tme(eset = exp,
                            method = "cibersort",
                            arrays = F,
                            perm = 1000
)
cell_bar_plot(input = cibersort[1:12,], features = colnames(cibersort)[3:24], title = "CIBERSORT Cell Fraction")

#cibersort<-cibersort[,-c(24:26)]

# IPS
ips <- deconvo_tme(eset = exp,
                      method = "ips",
                      plot = F
)
# quanTIseq
quantiseq <- deconvo_tme(eset = exp,
                         arrays = FALSE,
                         tumor = TRUE,
                            method = "quantiseq",
                            scale_mrna = T
)

# ESTIMATE

estimate <- deconvo_tme(eset = exp,
                           method = "estimate"
)
# TIMER
timer <- deconvo_tme(eset = exp
                        ,method = "timer"
                        ,group_list = rep("coad",dim(exp)[2])
)
# ssGSEA
load(file = "ssGSEA28.Rdata")
ssgsea <- calculate_sig_score(eset = exp
                                 , signature = cellMarker 
                                 , method = "ssgsea"
)

