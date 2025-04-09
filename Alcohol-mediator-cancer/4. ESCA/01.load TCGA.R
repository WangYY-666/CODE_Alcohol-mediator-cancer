
rm(list=ls())

library(dplyr)
library("rjson")

# 指定文件名和路径
metafile <- "ESCA/metadata.cart.2024-12-04.json"  # 下载的metadata文件名称
path1<-".\\ESCA\\Data\\"
outfilename <- "TCGA-ESCA.txt"  # 输出表达矩阵文件的名称

# 读取metadata文件
json <- jsonlite::fromJSON(metafile)
sample_id <- sapply(json$associated_entities, function(x) { x[, 1] })

# 创建包含样本ID和文件名的数据框
file_sample <- data.frame(sample_id, file_name = json$file_name)

# 列出所有包含基因计数的文件
count_file <-list.files(path1, pattern = 'gene_counts\\.tsv$', recursive = TRUE, full.names = FALSE)
# 提取文件名
count_file_name <- strsplit(count_file, split = '/')
count_file_name <- sapply(count_file_name, function(x) { x[2] })

# 循环读取每个基因计数文件，并将其合并到矩阵中
matrix = data.frame(matrix(nrow=60660,ncol=0)) 
matrix_tpm = data.frame(matrix(nrow=60660,ncol=0))

for (i in 1:length(count_file_name)){
  path <- paste0(path1, count_file[i])  # 构建文件路径
  data <- read.delim(path, fill = TRUE, header = FALSE, row.names = 1) 
  colnames(data)<-data[2,]
  data <-data[-c(1:6),]
  #取出COUNT矩阵:3，TPM:6，fpkm-unstranded:7，fpkm-up-unstranded:8
  data_tpm<-data[6]
  data <- data[3]
  f=unlist(strsplit(count_file_name[i],split='/')) 
  fn=f[length(f)]
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==fn)]
  colnames(data_tpm) <- file_sample$sample_id[which(file_sample$file_name==fn)]
  print(fn)
  matrix <- cbind(matrix,data)
  matrix_tpm<-cbind(matrix_tpm,data_tpm)
}

# 读取第一个样本的文件，用于获取基因名
sample1 <- paste0(path1, count_file[1])
names <- read.delim(sample1, fill = TRUE, header = FALSE, row.names = 1)
colnames(names) <- names[2, ]  # 设置列名为第二行内容
names <- names[-c(1:6), ]  # 删除前6行
names <- names[, 1:2]  # 选择前两列

# 获取矩阵和基因名的交集
same <- intersect(rownames(matrix), rownames(names))
same_tpm <- intersect(rownames(matrix_tpm), rownames(names))
# 筛选交集部分的数据
matrix <- matrix[same, ]
names <- names[same, ]

matrix_tpm <- matrix_tpm[same_tpm, ]
names_tpm <- names[same_tpm, ]
# 添加基因symbol到矩阵
matrix$symbol <- names[, 1]
matrix <- matrix[, c(ncol(matrix), 1:(ncol(matrix) - 1))]  # 调整列顺序

matrix_tpm$symbol <- names_tpm[, 1]
matrix_tpm <- matrix_tpm[, c(ncol(matrix_tpm), 1:(ncol(matrix_tpm) - 1))] 

# 输出结果到指定文件
write.table(matrix, file = outfilename, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(matrix_tpm, file = "TCGA-ESCA_TPM", row.names = FALSE, quote = FALSE, sep = "\t")

#将gene_name列去除重复的基因，保留每个基因最大表达量结果(较慢)
gene_name <-names[,1]
matrix0 <- aggregate( . ~ gene_name,data=matrix, max)   
exp <- matrix0[,-1]
rownames(exp) <- matrix0[,1]

#只保留表达量大于0的行并对每一列应用as.numeric函数
exp <- mutate_all(exp, function(x) as.numeric(as.character(x))) 
expr_ESCA <- exp[rowSums(exp)>0,]

save(expr_ESCA,file = "count_ESCA(clear).Rdata")




# 加载必要的库
library(readr)
library(tidyverse)
library(dplyr)
library(data.table)


# 读取clinical数据
clin <- fread("ESCA/clinical.tsv", data.table = FALSE)

# 查看数据的前几行，确保成功读取
head(clin)

# 提取感兴趣的列并删除重复项
clin_time <- clin %>%
  dplyr::select(
    case_submitter_id,
    vital_status,
    days_to_death,
    days_to_last_follow_up,
    age_at_index,
    gender,
    ajcc_pathologic_t,
    ajcc_pathologic_m,
    ajcc_pathologic_n,
    ajcc_pathologic_stage
  ) %>%
  dplyr::filter(!duplicated(case_submitter_id))


clin_merge <- clin_time %>%
  dplyr::mutate(
    OS.time = case_when(
      vital_status == "Alive" ~ days_to_last_follow_up,  # 若存活，使用 days_to_last_follow_up
      vital_status == "Dead" ~ days_to_death             # 若死亡，使用 days_to_death
    ),
    OS = case_when(
      vital_status == "Alive" ~ 0,  # 存活用0表示
      vital_status == "Dead" ~ 1    # 死亡用1表示
    )
  )


clin_merge1 <- clin_merge %>%
  dplyr::rename(
    age = age_at_index,       # 将age_at_index列名修改为age
    T = ajcc_pathologic_t,    # 将ajcc_pathologic_t列名修改为T
    M = ajcc_pathologic_m,    # 将ajcc_pathologic_m列名修改为M
    N = ajcc_pathologic_n,    # 将ajcc_pathologic_n列名修改为N
    stage = ajcc_pathologic_stage  # 将ajcc_pathologic_stage列名修改为stage
  )

# 查看数据的前几行以确保列名已修改
head(clin_merge1)

# 删除第3和第4列，并去除缺失值
clin <- clin_merge1[,-(3:4)]  # 删除第3和第4列
clin <- na.omit(clin)         # 删除包含NA的行

alcohol<- fread("ESCA/exposure.tsv", data.table = FALSE)
alcohol<-alcohol[,c(2,9)]
alcohol<-alcohol[alcohol$alcohol_history !="Not Reported",]
clin_ESCA<-merge(clin,alcohol,by="case_submitter_id")
colnames(clin_ESCA)[1]<-"id"
group<-condition_table[!duplicated(condition_table$id),]
clin_ESCA<-merge(group,clin_ESCA,by="id")
# 将结果保存为一个制表符分隔的文件（TSV格式）
write.table(clin_ESCA, file = "clinical_ESCA.txt", sep = "\t", row.names = FALSE, quote = FALSE)
save(clin_ESCA,file = "Clinic data(ESCA). Rdata")


###
load(file ="Clinic data(ESCA). Rdata")
clinc_TCGA<-data.table::fread(file = "ESCA/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix")
load(file = "count_ESCA(clear).Rdata")
clinc_need<-clinc_TCGA[,c(4,13,14,38,39,41,44,46:48,70,71,97,98)]
colnames(clinc_need)[1]<-"id"
clinc_need<-unique(clinc_need)
clinc_ESCA<-merge(clin_ESCA,clinc_need,by="id")
clinc_ESCA<-clinc_ESCA[,-14]

save(clinc_ESCA,file = "clinic data merged(ESCA).Rdata")


#___________________________________________________________________________________________#
###备用###
#1 载入临床特征数据，并进行其中缺失数据的清洗
clinical <- jsonlite::fromJSON('ESCA/clinical.cart.2024-12-04.json')
demographic_data <- as.data.frame(clinical$demographic)
demographic_data <- demographic_data[,c(3,5,7,9)]
colnames(demographic_data)[1:4] <- c("sex","fustat","age","futime")
clinical <- clinical[,4]
clinical <- cbind(clinical, demographic_data)
colnames(clinical)[1]<- "id" 

#2 区分肿瘤和正常样本
condition_table <- tidyr::separate(data=data.frame('ids'=colnames(exprSet)[1:length(exprSet)]),
                                   col='ids',
                                   sep='-',
                                   into=c('c1','c2','c3','c4','c5','c6','c7'))

condition_table <- cbind('ids'=colnames(exprSet)[1:length(exprSet)], condition_table) %>%
  cbind('submitter_IDs'=substr(colnames(exprSet)[1:length(exprSet)], 1, 12))
condition_table$c4 <- gsub(condition_table$c4, pattern = '^0..', replacement = 'cancer')
condition_table$c4 <- gsub(condition_table$c4, pattern = '^1..', replacement = 'normal')
condition_table <- condition_table[,c('ids','c4','submitter_IDs')]
names(condition_table) <- c('sample', 'group', 'submitter_id')
condition_table$submitter_id <- as.character(condition_table$submitter_id)
condition_table$sample <- as.character(condition_table$sample)
save(condition_table,file = 'condition_table')
colnames(condition_table)[3]<- "id" 
clinical <- condition_table %>%left_join(clinical, by = "id")



#3 输出最终文件
save(clinical,exprSet,file="TCGA初始数据count.RData")



