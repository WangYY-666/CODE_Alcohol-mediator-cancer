
rm(list = ls())

#1.获取数据
library(tinyarray)  
gse = "GSE93601"
geo= geo_download(gse)

library(limma)
#2.分组信息
pd<-geo$pd
exp<-geo$exp
exp=normalizeBetweenArrays(exp)
group_list = ifelse(pd$`tissue:ch1`=="breast tumor-adjacent normal","normal","tumor")
table(group_list)
group_list = factor(group_list,levels = c("normal","tumor")) 
#3.探针信息
# find_anno(gpl_number) 
gpl_data<-data.table::fread("BRCA/GPL22920_FFPE_annotations_NCBIupdated_18427probes_18417genesymbol_7Dec15.txt")
ids<-gpl_data[,c(1,3)]
colnames(ids)<-c("probe_id", "symbol" )
ids = ids[!duplicated(ids$symbol),]   #去除ids中重复的探针
ids<-as.data.frame(ids)

dcp = get_deg_all(exp,group_list,ids,symmetry = T,logFC_cutoff = 1 ,cluster_cols = F,entriz = F)

load(file = "BRCA data.Rdata")
clin_BRCA<-cbind(ID=rownames(clinc_data),clinc_data)
clin_BRCA1<-clin_BRCA[clin_BRCA$tissue=="breast tumor",]
exp_BRCA<-exp[,colnames(exp) %in% clin_BRCA1$ID]
p = identical(rownames(clin_BRCA1),colnames(exp_BRCA))
p
group= ifelse(clin_BRCA1$alcohol_history=="yes","alcohol_yes","alcohol_no")
table(group)
group= factor(group,levels = c("alcohol_yes","alcohol_no")) 
deg_alcohol<-get_deg(exp_BRCA,group,ids,
              logFC_cutoff=1,
              pvalue_cutoff=0.05,
              adjust = FALSE,
              entriz = FALSE,
              species = "human")
deg_alcohol_p<-deg_alcohol[deg_alcohol$P.Value<0.05,]
deg_alcohol_final<-deg_alcohol_p[order(abs(deg_alcohol_p$logFC), decreasing = TRUE)[1:10],]
deg_alcohol_final$change<-ifelse(deg_alcohol_final$logFC>0,"up","down")


x11()
dcp$plots
deg<-dcp$deg
rownames(deg)<-deg$symbol
deg<-deg[,-c(7,8)]
write.csv(deg,file = "01.normal-tumor/deg.csv")

deg1<-get_deg(exp_BRCA,group,ids,
                    logFC_cutoff=1,
                    pvalue_cutoff=0.05,
                    adjust = FALSE,
                    entriz = TRUE,
                    species = "human")

clin1<-clin_BRCA[clin_BRCA$alcohol_intake=="0 g/day",]
exp_1<-exp[,colnames(exp) %in% clin1$ID]
p = identical(rownames(clin1),colnames(exp_1))
p
group= ifelse(clin1$tissue=="breast tumor-adjacent normal","normal","tumor")
table(group)
group= factor(group,levels = c("normal","tumor")) 
deg2<-get_deg(exp_1,group,ids,
              logFC_cutoff=1,
              pvalue_cutoff=0.05,
              adjust = FALSE,
              entriz = TRUE,
              species = "human")
table(deg2$change)
clin2<-clin_BRCA[clin_BRCA$alcohol_intake=="0-10 g/day",]
exp_2<-exp[,colnames(exp) %in% clin2$ID]
p = identical(rownames(clin2),colnames(exp_2))
p
group= ifelse(clin2$tissue=="breast tumor-adjacent normal","normal","tumor")
table(group)
group= factor(group,levels = c("normal","tumor")) 
deg3<-get_deg(exp_2,group,ids,
              logFC_cutoff=1,
              pvalue_cutoff=0.05,
              adjust = FALSE,
              entriz = TRUE,
              species = "human")
table(deg3$change)
clin3<-clin_BRCA[clin_BRCA$alcohol_intake==">10 g/day",]
exp_3<-exp[,colnames(exp) %in% clin3$ID]
p = identical(rownames(clin3),colnames(exp_3))
p
group= ifelse(clin3$tissue=="breast tumor-adjacent normal","normal","tumor")
table(group)
group= factor(group,levels = c("normal","tumor")) 
deg4<-get_deg(exp_3,group,ids,
              logFC_cutoff=1,
              pvalue_cutoff=0.05,
              adjust = FALSE,
              entriz = TRUE,
              species = "human")
table(deg4$change)
save(deg1,deg2,deg3,deg4,file = "res_diffgene.Rdata")

library(clusterProfiler)
library(ggthemes)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)

# 1.GO 富集分析----
#deg<-dcp$deg
deg<-deg_alcohol_final
genes<-as.vector(deg$symbol)
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)       
entrezIDs <- as.character(entrezIDs)
deg=cbind(deg,ENTREZID=entrezIDs)

#(1)输入数据
gene_up = deg$ENTREZID[deg$change == 'up'] 
gene_down = deg$ENTREZID[deg$change == 'down'] 
gene_diff = c(gene_up,gene_down)
gene_diff<-gene_diff[gene_diff !="NA"]
#(2)富集
#以下步骤耗时很长，设置了存在即跳过
if(!file.exists(paste0(gse,"_GO.Rdata"))){
  ego <- enrichGO(gene = gene_diff,
                  OrgDb=org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "ALL",
                  #pAdjustMethod = "BH",
                  #minGSSize = 1,
                  #pvalueCutoff = 0.05,
                  #qvalueCutoff = 0.05,
                  readable = TRUE)
  #ont参数：One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
  save(ego,ego_BP,file = paste0(gse,"_GO.Rdata"))
}
ego_ALL<-as.data.frame(ego)
ego_ALL$Description
write.csv(ego_ALL,file = "./01.normal-tumor/GO_result.csv")
#(3)可视化
#条带图
p1<-barplot(ego)
library(tidyverse)
ggsave("./01.normal-tumor/GO条带图.pdf",p1)
#气泡图
p2<-dotplot(ego)
ggsave("./01.normal-tumor/GO气泡图.pdf",p2)

p3<-dotplot(ego, split = "ONTOLOGY", font.size = 10, 
            showCategory = 6) + 
  facet_grid(ONTOLOGY ~ ., space = "free_y",scales = "free_y") + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45))
ggsave("./01.normal-tumor/气泡拼接图.pdf",width = 20, height = 20, units = "cm")

#geneList 用于设置下面图的颜色
geneList = deg$logFC
names(geneList)=deg$ENTREZID

#(3)展示top通路的共同基因，要放大看。
#Gene-Concept Network
p4<-cnetplot(ego,categorySize="pvalue", foldChange=geneList,colorEdge = TRUE)
ggsave("./01.normal-tumor/top_pathway_gene01.pdf",width = 20, height = 20, units = "cm")
p5<-cnetplot(ego, showCategory = 3,foldChange=geneList, circular = TRUE, colorEdge = TRUE)
ggsave("./01.normal-tumor/top_pathway_gene02.pdf",width = 20, height = 20, units = "cm")
#Enrichment Map,这个函数最近更新过，版本不同代码会不同

#Biobase::package.version("enrichplot")

if(T){
  emapplot(pairwise_termsim(ego)) #新版本
}else{
  emapplot(ego)#旧版本
}
#(4)展示通路关系 https://zhuanlan.zhihu.com/p/99789859
#goplot(ego)
goplot(ego_BP)

#(5)Heatmap-like functional classification
heatplot(ego,foldChange = geneList,showCategory = 8)
??heatplot

##另一种画法
library(DOSE)
data(geneList)
de <- names(geneList)[1:50]
x <- enrichDO(de)
heatplot(x)

# 2.KEGG pathway analysis----
#上调、下调、差异、所有基因
#（1）输入数据
gene_up = deg[deg$change == 'up','ENTREZID'] 
gene_down = deg[deg$change == 'down','ENTREZID'] 
gene_diff = c(gene_up,gene_down)
gene_diff<-gene_diff[gene_diff !="NA"]
#（2）对上调/下调/所有差异基因进行富集分析
#install.packages('R.utils')
R.utils::setOption( "clusterProfiler.download.method",'auto' )
##解决办法见技能树公众号：KEGG数据库倒闭了吗

if(!file.exists(paste0(gse_number,"_KEGG.Rdata"))){
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa')
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa')
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa')
  #save(kk.diff,kk.down,kk.up,file = paste0(gse_number,"_KEGG.Rdata"))
}
load(paste0(gse_number,"_KEGG.Rdata"))

#(3)看看富集到了吗？https://mp.weixin.qq.com/s/NglawJgVgrMJ0QfD-YRBQg
table(kk.diff@result$p.adjust<0.05)
table(kk.up@result$p.adjust<0.05)
table(kk.down@result$p.adjust<0.05)
res_kegg<-as.data.frame(kk.diff@result)
write.csv(res_kegg,file = "./01.normal-tumor/KEGG_result.csv")
#(4)双向图
# 富集分析所有图表默认都是用p.adjust,富集不到可以退而求其次用p值，在文中说明即可
down_kegg <- kk.down@result %>%
  filter(pvalue<0.05) %>% #筛选行
  mutate(group=-1) #新增列

up_kegg <- kk.up@result %>%
  filter(pvalue<0.05) %>%
  mutate(group=1)

#source("kegg_plot_function.R")
kegg_plot <- function(up_kegg,down_kegg){
  dat=rbind(up_kegg,down_kegg)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  dat=dat[order(dat$pvalue,decreasing = F),]
  
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="-log10Pvalue") +
    coord_flip() + 
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    ggtitle("Pathway Enrichment") 
}
g_kegg <- kegg_plot(up_kegg,down_kegg)
g_kegg
#g_kegg +scale_y_continuous(labels = c(2,0,2,4,6))
ggsave(g_kegg,filename = './01.normal-tumor/kegg_up_down.pdf')
