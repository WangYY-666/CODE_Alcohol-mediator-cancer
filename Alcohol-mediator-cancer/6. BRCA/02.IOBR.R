.libPaths("D:/PackagesR/R_LIBS")
rm(list = ls())


library(IOBR)

load(file = "BRCA data.Rdata")
clin_BRCA<-cbind(ID=rownames(clinc_data),clinc_data)
clin_BRCA1<-clin_BRCA[clin_BRCA$tissue=="breast tumor",]
exp<-exprSet1[,colnames(exprSet1) %in% clin_BRCA1$ID]
dim(exp)
tme_deconvolution_methods

# MCPcounter
mcpcounter <- deconvo_tme(eset = exp,
                             method = "mcpcounter"
)
epic <- deconvo_tme(eset = exp,
                       method = "epic",
                       arrays = F
)
#epic<-epic[,-8]
# xCell
xcell <- deconvo_tme(eset = exp,
                        method = "xcell",
                        arrays = F
)
# CIBERSORT
cibersort <- deconvo_tme(eset = exp,
                            method = "cibersort",
                            arrays = F,
                            perm = 1000
)
save(cibersort,file = "GSE93601_cibersort.Rdata")
load(file = "GSE93601_cibersort.Rdata")
cibersort<-im_cibersort

cibersort<-cibersort[,-c(24:26)]
# IPS
ips <- deconvo_tme(eset = exp,
                      method = "ips",
                      plot = F
)

# quanTIseq
quantiseq <- deconvo_tme(eset = exp,
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


