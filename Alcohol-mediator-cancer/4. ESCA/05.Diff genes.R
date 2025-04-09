
rm(list=ls())


load(file = "count_ESCA(clear).Rdata")
load(file = "clinic_data(types).Rdata")

clin_ESCA<-clin_CA
clin_ESCA<-clin_SCC
library(dplyr)

rt1<-as.data.frame(t(expr_ESCA))
rt2<-cbind(sample=rownames(rt1),rt1)
rt3<-merge(clin_ESCA,rt2,by="sample")
rownames(rt3)<-rt3[,1]
#rt<-as.data.frame(t(rt3[,-c(1:13)]))
rt<-as.data.frame(t(rt3[,-c(1:26)]))

zero_percentage <- rowMeans(rt == 0)

rt_ESCA <- rt[zero_percentage < 0.5, ]
dim(rt_ESCA)

#save(rt_ESCA,file = "count_ESCA(deep clean)_CA.Rdata")
save(rt_ESCA,file = "count_ESCA(deep clean)_SCC.Rdata")
#差异分析

library(edgeR)
##创建分组情况
group <- factor(clin_ESCA$alcohol_history,levels = c("Yes","No"))
table(group)
group_list <- group
design <- model.matrix(~0+group_list)
rownames(design) <- colnames(rt_ESCA)
colnames(design) <- levels(group_list)

DGElist <- DGEList(counts = rt_ESCA, group = group_list)

keep_gene <- rowSums(cpm(DGElist) > 1 ) >= 2

DGElist <- DGElist[keep_gene, keep.lib.sizes = FALSE ]

DGElist <- calcNormFactors( DGElist )

DGElist <- estimateGLMCommonDisp(DGElist, design)

DGElist <- estimateGLMTrendedDisp(DGElist, design)

DGElist <- estimateGLMTagwiseDisp(DGElist, design)
fit <- glmFit(DGElist, design)
results <- glmLRT(fit, contrast = c(-1, 1))
ESCA_nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
ESCA_nrDEG_edgeR <- as.data.frame(ESCA_nrDEG_edgeR)


library(data.table) 
library(dplyr) 

ESCA_edgeR_DEG <- ESCA_nrDEG_edgeR
logFC_cutoff <- with(ESCA_edgeR_DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )
#logFC_cutoff <- 1
k1 = (ESCA_edgeR_DEG$PValue < 0.05)&(ESCA_edgeR_DEG$logFC < -logFC_cutoff)
k2 = (ESCA_edgeR_DEG$PValue < 0.05)&(ESCA_edgeR_DEG$logFC > logFC_cutoff)
ESCA_edgeR_DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))

ESCA_edgeR_DEG_signif <- ESCA_edgeR_DEG  %>% 
  filter(change !="NOT")

save(ESCA_edgeR_DEG,file = "ESCA_edgeR_DEG_1.Rdata")
save(ESCA_edgeR_DEG_signif,file = "ESCA_edgeR_DEG_signif_1.Rdata")




library(edgeR)
library(limma)
library(DESeq2)

Test_Expr <- rt_ESCA

Test_Expr = Test_Expr[rowMeans(Test_Expr)>1,] 

sample_yes<-clin_ESCA[clin_ESCA$alcohol_history=="Yes",]$sample
sample_no<-clin_ESCA[clin_ESCA$alcohol_history=="No",]$sample
alcohol_yes <- Test_Expr[,sample_yes]
alcohol_no  <- Test_Expr[,sample_no]
exprSet_by_group <- cbind(alcohol_yes,alcohol_no)
group_list <- c(rep('alcohol_yes',ncol(alcohol_yes)),rep('alcohol_no',ncol(alcohol_no)))

group_list = factor(group_list)
design <- model.matrix(~0+group_list)
rownames(design) = colnames(exprSet_by_group)
colnames(design) <- levels(group_list)

exp<-exprSet_by_group

dge <- DGEList(counts=exp)

keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2 ## 过滤
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
dge <- calcNormFactors(dge)

v <- voom(dge,design, normalize="quantile")
fit <- lmFit(v, design)

constrasts = paste(rev(levels(group_list)),collapse = "-")
cont.matrix <- makeContrasts(contrasts=constrasts,levels = design) 
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)

DEG = topTable(fit2, coef=constrasts, n=Inf)
DEG = na.omit(DEG)

logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )
#logFC_cutoff <- 1
k1 = (DEG$P.Value < 0.05)&(DEG$logFC < -logFC_cutoff)
k2 = (DEG$P.Value < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG$change)
head(DEG)
ESCA_limma_DEG <- DEG

ESCA_limma_DEG_signif <- ESCA_limma_DEG  %>% 
  filter(change !="NOT")

save(ESCA_limma_DEG,file = "ESCA_limma_DEG_1.Rdata")
save(ESCA_limma_DEG_signif ,file = "ESCA_limma_DEG_signif_1.Rdata")



library(DESeq2)

#Test_Expr = avereps(Test_Expr[,-1],ID = Test_Expr$gene_name) 
Test_Expr<-rt_ESCA

Test_Expr = Test_Expr[rowMeans(Test_Expr)>1,] 

sample_yes<-clin_ESCA[clin_ESCA$alcohol_history=="Yes",]$sample
sample_no<-clin_ESCA[clin_ESCA$alcohol_history=="No",]$sample

alcohol_yes <- Test_Expr[,sample_yes]
alcohol_no  <- Test_Expr[,sample_no]
exprSet_by_group <- cbind(alcohol_yes,alcohol_no)
group_list <- c(rep('alcohol_yes',ncol(alcohol_yes)),rep('alcohol_no',ncol(alcohol_no)))
#分组
condition = factor(group_list)
coldata <- data.frame(row.names = colnames(exprSet_by_group), condition)
dds <- DESeqDataSetFromMatrix(countData = exprSet_by_group,
                              colData = coldata,
                              design = ~condition)
dds$condition<- relevel(dds$condition, ref = "alcohol_no") # 指定对照组

dds <- DESeq(dds)
allDEG2 <- as.data.frame(results(dds))

logFC_cutoff <- with(allDEG2,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) )
#logFC_cutoff <- 1
k1 = (allDEG2$pvalue < 0.05)&(allDEG2$log2FoldChange < -logFC_cutoff)
k2 = (allDEG2$pvalue < 0.05)&(allDEG2$log2FoldChange > logFC_cutoff)
allDEG2$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(allDEG2$change)
head(allDEG2)
ESCA_DESeq2_DEG<-allDEG2
ESCA_DESeq2_DEG_signif <- ESCA_DESeq2_DEG  %>% 
  filter(change !="NOT")

save(ESCA_DESeq2_DEG,file = "ESCA_DESeq2_DEG_CA.Rdata")
save(ESCA_DESeq2_DEG_signif ,file = "ESCA_DESeq2_DEG_signif_CA.Rdata")


DEG_limma_voom <- ESCA_limma_DEG_signif
DEG_DEseq2 <- ESCA_DESeq2_DEG_signif
DEG_edgeR <-ESCA_edgeR_DEG_signif


allg <- intersect(rownames(DEG_limma_voom),rownames(DEG_edgeR))#取交集
allg <- intersect(allg,rownames(DEG_DEseq2))
ALL_DEG <- cbind(DEG_limma_voom[allg,c(1,4)],
                 DEG_DEseq2[allg,c(2,5)],
                 DEG_edgeR[allg,c(1,4)]) 

colnames(ALL_DEG)
colnames(ALL_DEG) <- c('limma_log2FC','limma_pvalue',
                       'DEseq2_log2FC','DEseq2_pvalue',
                       'edgeR_log2FC','edgeR_pvalue')

write.csv(ALL_DEG,file = 'DEG_results_CA.csv')
write.csv(ALL_DEG,file = 'DEG_results_SCC.csv')

colnames(ALL_DEG)
print(cor(ALL_DEG[,c(1,3,5)])) 

DEseq2_deg <- rownames(ESCA_DESeq2_DEG_signif)
edgeR_deg <- rownames(ESCA_edgeR_DEG_signif)
limma_deg <- rownames(ESCA_limma_DEG_signif)

gene_name <- as.data.frame(allg)
save(gene_name,file = "diffgenes_CA.Rdata")
save(gene_name,file = "diffgenes_SCC.Rdata")

library(Vennerable)

mylist <- list(DEseq2=DEseq2_deg, edgeR=edgeR_deg, limma=limma_deg)
str(mylist)
Vennplot <- Venn(mylist)
Vennplot1 <- Vennplot[,c('DEseq2','edgeR')]
Vennplot2 <- Vennplot[,c('DEseq2','limma')]
Vennplot3 <- Vennplot[,c('limma','edgeR')]

#pdf(file = paste0('Vennplot_lg2FC',log2FC_cutoff,'.pdf'))
plot(Vennplot, doWeights = F)
plot(Vennplot, doWeights = T) 
plot(Vennplot1, doWeights = T)
plot(Vennplot2, doWeights = T)
plot(Vennplot3, doWeights = T)
#dev.off()

#install.packages("UpSetR")
library(UpSetR)
upset(fromList(mylist),
      order.by = "freq",
      point.size =5,
      line.size =2,
      mainbar.y.label = " Count of Intersection",
      sets.x.label = "Datasets Size",
      text.scale =  c(1.5,1.5,1.5,1.5,2,1.5))
