#https://github.com/l-magnificence/Mime

rm(list=ls())

# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

#depens<-c('GSEABase', 'GSVA', 'cancerclass', 'mixOmics', 'sparrow', 'sva' , 'ComplexHeatmap' )
#for(i in 1:length(depens)){
#  depen<-depens[i]
#  if (!requireNamespace(depen, quietly = TRUE))  BiocManager::install(depen,update = FALSE)
#}

#if (!requireNamespace("CoxBoost", quietly = TRUE))
#  devtools::install_github("binderh/CoxBoost")

#if (!requireNamespace("fastAdaboost", quietly = TRUE))
#  devtools::install_github("souravc83/fastAdaboost")

#if (!requireNamespace("Mime1", quietly = TRUE))
#  devtools::install_github("l-magnificence/Mime")
set.seed(123)
library(Mime1)
library(tibble)
load(file = "diffgenes_HNSC.Rdata")
load(file = "TCGA-HNSC_tpm.RData")
genelist<-gene_name$allg
range(exprSet)
exprSet<-log2(exprSet+1)
exp1<-as.data.frame(t(exprSet))
exp1<-exp1[,colnames(exp1) %in% genelist]
exp1<-rownames_to_column(exp1,var = "sample")
exp<-merge(clin_HNSC[,c(2,5,4)],exp1,by="sample")
colnames(exp)[1]<-"ID"

train<-sample(exp$ID, dim(exp)[1]*0.6)
df1<-exp[exp$ID %in% train,]
test<- setdiff(exp$ID, train)
df2<-exp[exp$ID %in% test,]
list_train_vali_Data<-list(Dataset1=df1,Dataset2=df2)
#save(list_train_vali_Data,file = "seed123_list.Rdata")
load(file = "seed123_list.Rdata")
a<-list_train_vali_Data[["Dataset1"]][,-c(1:3)]
range(a)
#load("./External data/genelist.Rdata")
#load("./External data/Example.cohort.Rdata")
list_train_vali_Data[["Dataset1"]][1:5,1:5]
#list_train_vali_Data<-na.omit(list_train_vali_Data)
res <- Mime1::ML.Dev.Prog.Sig(train_data = list_train_vali_Data$Dataset1,
                              list_train_vali_Data = list_train_vali_Data,
                              unicox.filter.for.candi = T,
                              unicox_p_cutoff = 0.05,
                              candidate_genes = genelist,
                              mode = 'all',#all', 'single', and 'double'
                             # single_ml = c("RSF"),
                              nodesize =5,
                              seed = 123 
)
#save(res,file = "ML_res.Rdata")
p1<-cindex_dis_all(res,validate_set = names(list_train_vali_Data)[-1],order =names(list_train_vali_Data),width = 0.35)
pdf("ML_C_index.pdf", width = 8, height = 15) 
p1
dev.off()
cindex_dis_select(res,
                  model="RSF",
                  order= names(list_train_vali_Data))
survplot <- vector("list",2) 
for (i in c(1:2)) {
  print(survplot[[i]]<-rs_sur(res, model_name = "StepCox[forward] + RSF",dataset = names(list_train_vali_Data)[i],
                              #color=c("blue","green"),
                              median.line = "hv",
                              cutoff = 0.5,
                              conf.int = T,
                              xlab="Day",pval.coord=c(1000,0.9)))
}
aplot::plot_list(gglist=survplot,ncol=2)
all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Dataset1"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 1,
                             auc_cal_method="KM")
all.auc.2y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Dataset1"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 2,
                             auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Dataset1"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 3,
                             auc_cal_method="KM")
p6<-auc_dis_all(all.auc.1y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=1)
p7<-auc_dis_all(all.auc.2y,
                dataset = names(list_train_vali_Data),
                validate_set=names(list_train_vali_Data)[-1],
                order= names(list_train_vali_Data),
                width = 0.35,
                year=2)
p8<-auc_dis_all(all.auc.3y,
                dataset = names(list_train_vali_Data),
                validate_set=names(list_train_vali_Data)[-1],
                order= names(list_train_vali_Data),
                width = 0.35,
                year=3)
pdf("ML_AUC_1.pdf", width = 8, height = 15) 
p6
dev.off()
pdf("ML_AUC_2.pdf", width = 8, height = 15) 
p7
dev.off()
pdf("ML_AUC_3.pdf", width = 8, height = 15) 
p8
dev.off()

auc_dis_select(list(all.auc.1y,all.auc.2y,all.auc.3y),
               model_name="RSF",
               dataset = names(list_train_vali_Data),
               order= names(list_train_vali_Data),
               year=c(1,2,3))
roc_vis(all.auc.1y,
        model_name = "RSF",
        dataset = names(list_train_vali_Data),
        order= names(list_train_vali_Data),
        anno_position=c(0.65,0.55),
        year=1)
