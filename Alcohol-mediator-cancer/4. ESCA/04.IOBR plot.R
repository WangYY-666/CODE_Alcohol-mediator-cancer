
rm(list=ls())

library(IOBR)
library(tidyverse)
library(clusterProfiler)
library(xCell)

load(file = "Clinic data(ESCA). Rdata")


rt_tpm<-data.table::fread(file = "TCGA-ESCA_TPM")
rt_tpm<-as.data.frame(rt_tpm)
rownames(rt_tpm)<-rt_tpm[,1]


gene_name <-rt_tpm$symbol
matrix0 <- aggregate( . ~ gene_name,data=rt_tpm, max)   
rt1 <- matrix0[,-c(1,2)]
rownames(rt1) <- matrix0[,1]

#只保留表达量大于0的行并对每一列应用as.numeric函数
rt1 <- mutate_all(rt1, function(x) as.numeric(as.character(x))) 
rt1 <- rt1[rowSums(rt1)>0,]
#识别异常样本
ddd<-find_outlier_samples(eset = rt1,  show_plot = TRUE)
rt2<-rt1[,!(names(rt1) %in% ddd)]


rt2<-as.data.frame(t(rt2))
rt2<-rt2[-1,]
rt2<-cbind(rownames(rt2),rt2)
colnames(rt2)[1]<-"sample"
rt3<-merge(clin_ESCA,rt2,by="sample")
rt4<-rt3[rt3$group =="cancer",]
#save(rt4,file = "ESCA_tpm_clin (clean).Rdata")
rt_tpm<-rt4
save(rt_tpm,file = "ESCA_tpm_clin (176 samples).Rdata")


load(file = "ESCA_tpm_clin (clean).Rdata")
exp<-rt_tpm[,-c(1:13)]
rownames(exp)<-rt_tpm[,1]
exp<-as.data.frame(t(exp))
exp_numeric <- data.frame(lapply(exp, function(x) {
  if (is.factor(x)) {
    as.numeric(as.character(x))  # 对因子类型进行特殊处理
  } else {
    suppressWarnings(as.numeric(as.character(x)))  # 对其他类型尝试转换，并抑制可能的警告
  }
}))
rownames(exp_numeric)<-rownames(exp)
colnames(exp_numeric)<-colnames(exp)

#exp_log<-log2(exp_numeric+0.1)
#exp_log<-exp_log[apply(exp_log,1,sd)>0.5,]



tme_deconvolution_methods
cibersort<-deconvo_tme(eset = exp_numeric, method = "cibersort", arrays = TRUE, perm = 100 )
cibersort0<-deconvo_tme(eset = exp, method = "cibersort", arrays = TRUE, perm = 100 )

# head(cibersort)
cell_bar_plot(input = cibersort[1:12,], features = colnames(cibersort)[2:23], title = "CIBERSORT Cell Fraction")

# help(deconvo_epic)
epic <- deconvo_tme(eset = exp_numeric,
                    method = "epic",
                    arrays = F
)
mcp<-deconvo_tme(eset = exp_numeric, method = "mcpcounter")
xcell<-xCellAnalysis(exp_numeric,parallel.sz=1)
xcell<-as.data.frame(t(xcell))
xcell<-cbind(ID=rownames(xcell),xcell)
estimate<-deconvo_tme(eset = exp_numeric, method = "estimate")
timer<-deconvo_tme(eset = exp_numeric, method = "timer", group_list = rep("ESCA",dim(exp_numeric)[2]))
quantiseq<-deconvo_tme(eset = exp_numeric, tumor = TRUE, arrays = TRUE, scale_mrna = TRUE, method = "quantiseq")
ips<-deconvo_tme(eset = exp_numeric, method = "ips", plot= FALSE)

tme_combine<-xcell %>% 
  inner_join(.,cibersort[,-c(24:26)],by     = "ID") %>%
  inner_join(.,quantiseq,by = "ID") %>% 
  inner_join(.,mcp,by       = "ID") %>% 
  inner_join(.,epic,by      = "ID") %>% 
#  inner_join(.,estimate,by  = "ID") %>% 
  inner_join(.,timer,by     = "ID")
 # inner_join(.,ips,by       = "ID")
dim(tme_combine)

clin_ESCA1<-clin_ESCA[,c(2,13)]
colnames(clin_ESCA1)<-c("ID","Alcohol_history")
tme<-merge(clin_ESCA1,tme_combine,by="ID")
tme_mid<-tme[,-1]
rownames(tme_mid)<-tme[,1]

result_tme_comb = as.data.frame(t(tme_mid))

# sample_group = as.data.frame(t(result_tme_comb[2,]))

heat_data = result_tme_comb[-1,]

gene_group<-data.frame(cell =rownames(heat_data))%>%
  mutate(cell_group = c(rep('xcell',67),rep('cibersort',22),
                        rep('quantiseq',11),rep('mcpcounter',10),
                        rep("epic",8),rep("timer",6)))

for(i in names(heat_data)){heat_data[,i]<-as.numeric(heat_data[,i])}

heat_data = heat_data%>%
  mutate(cell_group=gene_group$cell_group)%>%
  dplyr::select(cell_group,everything())

str(heat_data)
cli = clin_ESCA%>%
  dplyr::select(c(colnames(.)[c(2,5:12)],"alcohol_history"))%>%
  mutate(gender=as.factor(gender),Stage=as.factor(stage),alcohol_history=as.factor(alcohol_history),
         stage_T=as.factor(T),stage_N=as.factor(N),stage_M=as.factor(M),
         OS=factor(OS,levels=c('0','1'),labels=c('alive','deaed')))
cli$OS.time<-as.numeric(cli$OS.time)
str(cli)
result_tidy<-
  heat_data %>%
  #这里将行名命名为新的列
  as_tibble(rownames = 'cell') %>%
  #这里选择去掉哪些列来进行归一化
  mutate_at(vars(-colnames(gene_group)), scale) %>%
  #宽数据转长数据
  pivot_longer(cols = -colnames(gene_group), 
               names_to = "sample", 
               values_to = "Value")%>%
  left_join(cli,by='sample')

library(ggsci)
library(RColorBrewer)
library(ComplexHeatmap)
pal_anno1<-pal_nejm()(7)

pal_plot1<-
  #调色板指定范围
  #seq范围越小颜色越深
  circlize::colorRamp2(seq(-3, 3, length.out = 3), 
                       # rev(RColorBrewer::brewer.pal(9, "RdBu"))
                       c('steelblue','white','firebrick')
  )
result_tidy1 = result_tidy%>%#'_CIBERSORT|_xCell|_quantiseq|_MCPcounter'
  mutate(cell=gsub('_xCell','',cell))%>%
  mutate(cell=gsub('_CIBERSORT',' ',cell))%>%
  mutate(cell=gsub('_quantiseq','  ',cell))%>%
  mutate(cell=gsub('_MCPcounter','   ',cell))%>%
  mutate(cell=gsub('_EPIC','    ',cell))%>%
  mutate(cell=gsub('_TIMER','     ',cell))
plot_heatmap1<-
  result_tidy1 %>%
  #这里指定分类用的列，包括行分类和列分类
  group_by(`cell_group`,`alcohol_history`) %>%
  heatmap(
    .row=cell, 
    .column=sample, 
    .value=Value,
    .scale = "row",
    #这里指定分类用的颜色板
    palette_grouping = list(pal_anno1,#列的颜色注释
                            c("#66C2A5", "steelblue")),#行的颜色注释
    #这里指定画图用的颜色板
    palette_value = pal_plot1,
    # palette_value = c("green", 
    #                   "black", "red")
    cluster_rows = F,
    border=T,
    row_names_gp = gpar(fontsize = 8),#设置行名字体大小
    show_column_names = F,
    column_order=rownames(tme_mid),#设置列及样本的顺序
    
  )

plot_heatmap1

col1<-c('#fee5d9','#fcae91','#fb6a4a','#cb181d')
col2<-c('#f2f0f7','#cbc9e2','#9e9ac8','#6a51a3')
plot_heatmap = plot_heatmap1%>%
  add_tile(`Stage`,palette = col1)%>%
  add_tile(`stage_T`,palette = col2)%>%
  add_tile(`stage_N`,palette = col2)%>%
  add_tile(`stage_M`,palette = col2)%>%
  add_tile(`gender`,palette = c('pink','skyblue'))%>%
  add_tile(`OS`,palette=c('white','black'))%>%
  add_tile(`OS.time`)%>%
  add_point(`age`)

plot_heatmap

pdf('./Heatmap_ESCA1.pdf',width = 16,height = 16)
print(plot_heatmap)
dev.off()


library(ggplot2)
colnames(estimate)[1]<-"sample"
df<-merge(id[,c(1,13)],estimate,by="sample")
df<-merge(clin_SCC[,c(2,13)],estimate,by="sample")
# 创建一个函数来绘制单个评分的箱线图
plot_boxplot <- function(df, score_col) {
  ggplot(df, aes_string(x = "alcohol_history", y = score_col, fill = "alcohol_history")) +
    geom_boxplot() +
    labs(title = paste("Boxplot for", score_col),
         x = "alcohol_history",
         y = score_col) +
    theme_minimal() +
    scale_fill_manual(values = c("No" = "blue", "Yes" = "red")) # 自定义填充颜色
}

# 创建一个空的图形对象，用于后续添加子图
p <- ggplot() + theme_void()

# 使用gridExtra包来合并多个图形（如果未安装，请先安装）
# install.packages("gridExtra")
library(gridExtra)

# 创建一个图形列表，用于存储每个评分的箱线图
plots <- list()

# 为每个评分绘制箱线图，并添加到图形列表中
for (score in c("StromalScore_estimate", "ImmuneScore_estimate", 
                "ESTIMATEScore_estimate", "TumorPurity_estimate")) {
  plots[[score]] <- plot_boxplot(df, score)
}

# 使用grid.arrange函数将多个图形合并为一个图形
# 2行2列布局，根据评分数量调整布局参数
grid.arrange(grobs = plots, ncol = 2, nrow = 2)
