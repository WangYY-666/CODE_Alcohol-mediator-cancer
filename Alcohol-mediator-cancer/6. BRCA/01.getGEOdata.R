
rm(list=ls())


library(GEOquery)
gse_number = "GSE93601"
eSet <- getGEO(gse_number,
               destdir = '.', 
               getGPL = F)
class(eSet)
length(eSet)
eSet = eSet[[1]]

exp<- exprs(eSet)
dim(exp)

pd <- pData(eSet)

p = identical(rownames(pd),colnames(exp))
p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]


ex <- exp    
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) {
  ex[which(ex <= 0)] <- NaN
  exp <- log2(ex)
  print("log2 transform finished")}else{
    print("log2 transform not needed")}
library(limma) 
boxplot(exp,outline=FALSE, notch=T, las=2)
exp=normalizeBetweenArrays(exp)
boxplot(exp,outline=FALSE, notch=T, las=2)
exp <- as.data.frame(exp)


gpl_data<-data.table::fread("BRCA/GPL22920_FFPE_annotations_NCBIupdated_18427probes_18417genesymbol_7Dec15.txt")
#gpl_data<-data.table::fread("BRCA/GPL22920_hGlue3_0.r1.TC_Annov.txt")
exp1<-cbind(probeset_id=rownames(exp),exp)
exp1<-as.data.frame(exp1)
expr<-merge(gpl_data,exp1,by="probeset_id")

exprSet1<-expr[,-c(1,2,4)]
colnames(exprSet1)[1]<-"Gene Symbol"
gene <- exprSet1$`Gene Symbol`
symbol <- strsplit(gene,split = " /// ",fixed = T)
Gene <- sapply(symbol,function(x){x[1]})
exprSet1$`Gene Symbol` <- Gene
exprSet1 <- as.data.frame(exprSet1)


rownames(exprSet1) <- exprSet1$`Gene Symbol` 
exprSet1 <- distinct(exprSet1,`Gene Symbol`,.keep_all = T)
rownames(exprSet1) <- exprSet1$`Gene Symbol`  
exprSet1 <- na.omit(exprSet1)               
rownames(exprSet1) <- exprSet1$`Gene Symbol`  
exprSet1 <- exprSet1[,-1]


colnames(pd)
clinc_data<-pd[,c(10,45:49)]
clinc_data<-clinc_data[clinc_data$`alcohol intake:ch1`!="NA",]
table(clinc_data$`alcohol intake:ch1`)
library(dplyr)
clinc_data<-clinc_data %>%
  mutate(alcohol_history = case_when(
    `alcohol intake:ch1` == ">10 g/day" ~ "yes",
    `alcohol intake:ch1` == "0-10 g/day" ~ "yes",
    `alcohol intake:ch1` == "0 g/day" ~ "no"))
colnames(clinc_data)<-c("tissue","diagnosis_age","alcohol_intake","avg_bodysize_age","bmi","ER receptor", "alcohol_history")
table(clinc_data$alcohol_history)
table(clinc_data$tissue)
library(stringr)
clinc_data$tissue<-str_remove(clinc_data$tissue, "^tissue: ")


p = identical(rownames(clinc_data),colnames(exprSet1))
p
if(!p) exprSet1 = exprSet1[,match(rownames(clinc_data),colnames(exprSet1))]



save(exprSet1,clinc_data,file = "BRCA data.Rdata")
