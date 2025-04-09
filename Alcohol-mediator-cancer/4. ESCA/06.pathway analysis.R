
rm(list=ls())


library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

load(file = "res_p_Single_Cox_SCC.Rdata")
load(file = "res_p_Single_Cox_CA.Rdata")
#rt<-gene_name$allg
rt<-rownames(res_p)
genes<-as.vector(rt)
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)       
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
colnames(out)=c("symbol","entrezID")
out=out[is.na(out[,"entrezID"])==F,]
out<-as.data.frame(out)
gene=out$entrezID
gene<-gene[gene != "NA"]

ego_ALL <- enrichGO(gene = gene,
                    OrgDb=org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    minGSSize = 1,
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)

ego_CC <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_BP <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_MF <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)


ego_ALL <- as.data.frame(ego_ALL)
ego_result_BP <- as.data.frame(ego_BP)
ego_result_CC <- as.data.frame(ego_CC)
ego_result_MF <- as.data.frame(ego_MF)

ego_result_BP <- ego_result_BP[order(ego_result_BP$Count,decreasing = T),]
ego_result_CC <- ego_result_CC[order(ego_result_CC$Count,decreasing = T),]
ego_result_MF <- ego_result_MF[order(ego_result_MF$Count,decreasing = T),]
ego <- rbind(ego_result_BP,ego_result_CC,ego_result_MF)
ego_ALL$Description
ego_result_BP$Description
ego_result_CC$Description
ego_result_MF$Description


#SCC
write.csv(ego_ALL,file = "pathway/SCC/ego_ALL.csv",row.names = T)
write.csv(ego_result_BP,file = "pathway/SCC/ego_result_BP.csv",row.names = T)
write.csv(ego_result_CC,file = "pathway/SCC/ego_result_CC.csv",row.names = T)
write.csv(ego_result_MF,file = "pathway/SCC/ego_result_MF.csv",row.names = T)


go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),                         
  Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("biological process", nrow(ego_result_BP)), 
                rep("cellular component", nrow(ego_result_CC)),
                rep("molecular function", nrow(ego_result_MF))), 
              levels=c("biological process", "cellular component","molecular function" )))
for(i in 1:nrow(go_enrich_df)){
  description_splite=strsplit(go_enrich_df$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") 
  go_enrich_df$Description[i]=description_collapse
  go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
}


go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
library(RColorBrewer)
COLS <- brewer.pal(name="Dark2",3)

ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + 
  geom_bar(stat="identity", width=0.9) + 
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The GO Terms of ESCC with Alcohol History")+
  theme_bw()


library(tidyverse)
library(ggfun)
library(grid)

Down_result <- ego_ALL %>%
  dplyr::select(1,3,6,9,13) %>%
  dplyr::group_by(ONTOLOGY) %>%
  dplyr::slice(1,3,5,6,8:10,11,13:15,21:23) %>%
  #dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()
p <- Down_result %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered = T)) %>%
  ggplot(aes(x = Count, y = Description)) + 
  geom_bar(aes(fill = pvalue), stat = "identity", width = 0.75) + 
  geom_text(aes(x = Count + 2, label = str_c(Count, str_c("(",sprintf("%.2e", pvalue), ")", sep = ""), sep = " "))) + 
  ggtitle(label = "GO enrichment of ESCC with Alcohol History") + 
  labs(x = "Count", y = "Description") + 
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) + 
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = "#000000"),
    axis.title = element_text(size = 15, color = "#000000"),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 0.75),
    legend.background = element_roundrect(color = "#969696"),
    plot.title = element_text(size = 15, hjust = 0.5)
  ) + 
  coord_cartesian(clip = "off")

p

p + 
  annotation_custom(
    grob = rectGrob(
      x = unit(1, "mm"),
      y = unit(1, "mm"),
      width = unit(45, "mm"),
      height = unit(8, "mm"),
      gp = gpar(fill = "#6baed6", col = NA, alpha = 0.5)
    ),
    xmin = unit(-2, "native"),
    xmax = unit(-2, "native"),
    ymin = unit(0.9, "native"),
    ymax = unit(1.0, "native")
  )  + 
  annotation_custom(
    grob = rectGrob(
      x = unit(1, "mm"),
      y = unit(1, "mm"),
      width = unit(110, "mm"),
      height = unit(8, "mm"),
      gp = gpar(fill = "#6baed6", col = NA, alpha = 0.5)
    ),
    xmin = unit(-42, "native"),
    xmax = unit(-42, "native"),
    ymin = unit(14.9, "native"),
    ymax = unit(15.0, "native")
  )


library(limma)
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(ggplot2)
library(RColorBrewer)


diff <- genes

columns(org.Hs.eg.db)

diff_entrez <- bitr(diff,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = "org.Hs.eg.db")
head(diff_entrez) 
KEGG_diff <- enrichKEGG(gene = diff_entrez$ENTREZID,
                        organism = "hsa", 
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        minGSSize = 10,
                        maxGSSize = 500)


KEGG_diff <- setReadable(KEGG_diff,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID")

KEGG_diff2 <- mutate(KEGG_diff,
                     RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))


KEGG_diff2 <- mutate(KEGG_diff2, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
KEGG_diff2@result$RichFactor[1:6]
KEGG_diff2@result$FoldEnrichment[1:6]

KEGG_result <- KEGG_diff2@result
KEGG_result$Description

