
rm(list=ls())

load(file = "BRCA data.Rdata")
clinc_data<-cbind(ID=rownames(clinc_data),clinc_data)
library(limma)
phe=clinc_data

table(phe$source_name_ch1)

group_list = ifelse(clinc_data$tissue=="breast tumor-adjacent normal","normal","tumor")
table(group_list)
dat[1:4,1:4]
dat<-as.data.frame(exprSet1)
dat[1:4,1:4]
dim(dat)
dat['ACTB',1]
dat['GAPDH',1]
#View(dat)
dat[1,]
dat[2,]

library(patchwork)
library(cowplot)
options(stringsAsFactors = F)
table(group_list)

dat[1:4,1:4]

dat=t(dat)
dat=as.data.frame(dat)
library("FactoMineR")
library("factoextra")
dat.pca <- PCA(dat , graph = FALSE)
p1<-fviz_pca_ind(dat.pca,
                 geom.ind = "point", 
                 col.ind =  group_list, 
                 addEllipses = T,
                 legend.title = "Groups"
)
p1

cg=names(tail(sort(apply(dat,1,sd)),1000))
library(pheatmap)
pheatmap(dat[cg,],show_colnames =F,show_rownames = F) 

n=scale(t(dat[cg,])) 
n[n>2]=2
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)

ac=data.frame(group=group_list)
n<-as.data.frame(t(n))
rownames(ac)=colnames(n)


pheatmap(n,show_colnames =F,show_rownames = F,
         annotation_col=ac)
p2<-pheatmap(n,show_colnames =F,show_rownames = F,
             annotation_col=ac,filename = 'heatmap_top1000_sd.png')

require(ggplotify)
p2 = as.ggplot(p2)
p1+p2
save(p1,p2,file = "patchwork.Rdata")

ggsave(p2,filename = './outputs/heatmap_top1000_sd.png',height = 5,width = 5)
dev.off()

dat[1:4,1:4] 
exprSet=dat
pheatmap::pheatmap(cor(exprSet)) 

colD=data.frame(group=group_list)
rownames(colD)=colnames(exprSet)
pheatmap::pheatmap(cor(exprSet),
                   annotation_col = colD,
                   show_rownames = F)
p3<-pheatmap::pheatmap(cor(exprSet),
                       annotation_col = colD,
                       show_rownames = F,
                       filename = 'cor_all.png')
ggsave(p3,filename = './outputs/heatmap_cor_all.png',height = 5,width = 5)
dev.off()

exprSet=exprSet[names(sort(apply(exprSet, 1,mad),decreasing = T)[1:10]),]
dim(exprSet) 
M=cor(exprSet)
pheatmap::pheatmap(M,annotation_col = colD)
pheatmap::pheatmap(M,
                   show_rownames = F,
                   annotation_col = colD)
p4<-pheatmap::pheatmap(M,
                       show_rownames = F,
                       annotation_col = colD,
                       filename = 'cor_top10.png')
ggsave(p4,filename = './outputs/cor_top10.png',height = 5,width = 5)
dev.off()

dat[1:4,1:4] 
table(group_list) 

boxplot(dat[1,]~group_list) 

library(limma)
design=model.matrix(~factor( group_list ))
design
fit=lmFit(dat,design)
fit=eBayes(fit)

options(digits = 4) 
#topTable(fit,coef=2,adjust='BH') 
deg = topTable(fit,coef=2,adjust='BH', n=Inf) 

head(deg) 

library(dplyr)
deg <- mutate(deg,symbol=rownames(deg))
head(deg)


library(dplyr)
deg_up <- deg %>% filter(logFC>1,P.Value<0.05)
deg_down <- deg %>% filter(logFC<=-1,P.Value<0.05)

write.csv(head(deg_up,10),file = "outputs/up10gene.csv")
write.csv(head(deg_down,10),file = "outputs/down10gene.csv")


logFC_t=1
P.Value_t = 0.05
k1 = (deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value < P.Value_t)&(deg$logFC > logFC_t)
deg <- mutate(deg,change = ifelse(k1,"down",ifelse(k2,"up","stable")))
table(deg$change)

library(clusterProfiler)
library(org.Hs.eg.db)
head(deg)
s2e <- bitr(deg$symbol, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)#人类

head(s2e)
dim(deg)
deg <- inner_join(deg,s2e,by=c("symbol"="SYMBOL"))
dim(deg)
length(unique(deg$symbol))
save(group_list,deg,logFC_t,P.Value_t,gse_number,file = "step3output.Rdata")

library(EnhancedVolcano)
p3<-EnhancedVolcano(deg,
                    lab =  rownames(deg),
                    x = 'logFC',
                    y = 'P.Value')
ggsave(filename = './outputs/EnhancedVolcano_for_DEG_DEseq2.pdf',width = 15,height = 15)
load("patchwork.Rdata")
p1+p2+p3
ggsave("./outputs/patchwork.png",width = 40, height = 30, units = "cm")

bp(dat[rownames(deg)[1],])
## for volcano 
if(T){
  nrDEG=deg
  head(nrDEG)
  attach(nrDEG)
  plot(logFC,-log10(P.Value))
  library(ggpubr)
  df=nrDEG
  df$v= -log10(P.Value) #df新增加一列'v',值为-log10(P.Value)
  ggscatter(df, x = "logFC", y = "v",size=0.5)
  
  df$g=ifelse(df$P.Value>0.01,'stable', #if 判断：如果这一基因的P.Value>0.01，则为stable基因
              ifelse( df$logFC >2,'up', #接上句else 否则：接下来开始判断那些P.Value<0.01的基因，再if 判断：如果logFC >1.5,则为up（上调）基因
                      ifelse( df$logFC < -2,'down','stable') )#接上句else 否则：接下来开始判断那些logFC <1.5 的基因，再if 判断：如果logFC <1.5，则为down（下调）基因，否则为stable基因
  )
  table(df$g)
  df$name=rownames(df)
  head(df)
  ggscatter(df, x = "logFC", y = "v",size=0.5,color = 'g')
  p4<-ggscatter(df, x = "logFC", y = "v", color = "g",size = 0.5,
                label = "name", repel = T,
                #label.select = rownames(df)[df$g != 'stable'] ,
                label.select = c('TTC9', 'AQP3', 'CXCL11','PTGS2'), #挑选一些基因在图中显示出来
                palette = c("#00AFBB", "#E7B800", "#FC4E07") )
  
  ggsave(p4,filename = './outputs/volcano.png',width = 15,height = 15)
  
  ggscatter(df, x = "AveExpr", y = "logFC",size = 0.2)
  df$p_c = ifelse(df$P.Value<0.001,'p<0.001',
                  ifelse(df$P.Value<0.01,'0.001<p<0.01','p>0.01'))
  table(df$p_c )
  p5<-ggscatter(df,x = "AveExpr", y = "logFC", color = "p_c",size=0.2, 
                palette = c("green", "red", "black") )
  
  ggsave(p5,filename = './outputs/MA.png',width = 15,height = 15)  
  
}

rm(list=ls())
## for heatmap 
if(T){ 
  load("step3output.Rdata")
  load(file = 'step1-output.Rdata')
  dat[1:4,1:4]
  table(group_list)
  x=deg$logFC 
  head(deg)
  names(x)=deg$symbol 
  cg=c(names(head(sort(x),5)),
       names(tail(sort(x),5)))
  library(pheatmap)
  ac=data.frame(group=group_list)
  rownames(ac)=colnames(dat[cg,])
  p<-pheatmap(dat[cg,],show_colnames =T,show_rownames = T,
              color = colorRampPalette(colors = c("blue","white","red"))(100),
              annotation_legend = T,annotation_col=ac)
  ggsave('./outputs/heatmap_TOP10_hub_genes.png',p)}
