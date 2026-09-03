library(readr)      
library(vroom)      
library(tidyr)      
library(tibble)     
library(dplyr)      
library(TwoSampleMR) 
library(ggplot2)   
library(gt)         
library(meta)     
library(ComplexHeatmap) 
library(reshape2)   
library(grid)       
library(dplyr)

res_group<-read.table(file="alcohol drinkers status/alcohol group.txt",sep="\t",header = TRUE)
colnames(res_group)[2]<-"id.exposure"
res_ivw<-read.csv(file = "alcohol drinkers status/res_ivw.csv")

cancer_list <- list(
  BLCA = list(name = "Bladder cancer", ids = c("ieu-b-4874")),
  BRCA = list(name = "Breast cancer", ids = c("ebi-a-GCST90018799", "ieu-a-1168", "ieu-b-4810", "ukb-a-55", "ukb-b-16890")),
  CESC = list(name = "Cervical cancer", ids = c("ebi-a-GCST90018817")),
  CHOL = list(name = "Hepatic bile duct cancer", ids = c("ebi-a-GCST90018803")),
  COAD = list(name = "Colon adenocarcinoma", ids = c("finn-b-C3_COLON_ADENO_EXALLC")),
  DLBC = list(name = "Diffuse large B-cell lymphoma", ids = c("finn-b-C3_DLBCL_EXALLC")),
  ESCA = list(name = "Esophageal cancer", ids = c("ebi-a-GCST90018841","ebi-a-GCST90018621","bbj-a-117")),
  GBM = list(name = "Brain glioblastoma", ids = c("finn-b-C3_GBM_EXALLC")),
  HNSC = list(name = "Head and neck cancer", ids = c("ieu-b-4912")),
  KIPAN = list(name = "Kidney cancer", ids = c("finn-b-C3_KIDNEY_NOTRENALPELVIS_EXALLC", "finn-b-CD2_BENIGN_KIDNEY_EXALLC")),
  LAML = list(name = "Acute myeloid leukaemia", ids = c("finn-b-C3_AML_EXALLC")),
  LIHC = list(name = "Liver cell carcinoma", ids = c("ieu-b-4953")),
  LUAD = list(name = "Lung adenocarcinoma", ids = c("ieu-a-965")),
  LUSC = list(name = "Squamous cell lung carcinoma", ids = c("ebi-a-GCST004750")),
  MESO = list(name = "Mesothelioma", ids = c("finn-b-C3_MESOTHELIOMA_EXALLC")),
  OV = list(name = "Ovarian cancer", ids = c("ieu-b-4963", "ieu-a-1120","ieu-a-1237")),
  PAAD = list(name = "Pancreatic cancer", ids = c("ebi-a-GCST90018893","bbj-a-140")),
  PRAD = list(name = "Prostate cancer", ids = c("ieu-b-4809")),
  READ = list(name = "Rectal cancer", ids = c("finn-b-C3_RECTUM_EXALLC", "finn-b-CD2_BENIGN_RECTUM_EXALLC")),
  SKCM = list(name = "Melanoma skin cancer", ids = c("ieu-b-4969")),
  STAD = list(name = "Stomach cancer", ids = c("finn-b-CD2_BENIGN_STOMACH_EXALLC", "finn-b-C3_STOMACH_EXALLC")),
  TGCT = list(name = "Testicular cancer", ids = c("finn-b-C3_TESTIS_EXALLC")),
  THCA = list(name = "Thyroid cancer", ids = c("ebi-a-GCST90018929")),
  UCEC = list(name = "Endometrial cancer", ids = c("ukb-b-13545"))
)

#dir.create("02. meta", showWarnings = FALSE)

meta_results <- data.frame()


TCGA_CANCER <- c("BLCA","BRCA","CESC", "CHOL","COAD","DLBC",'ESCA',
                 "GBM","HNSC","KIPAN","LAML","LIHC", "LUAD","LUSC",
                 "MESO","OV","PAAD","PRAD","READ", "SKCM","STAD",
                 "TGCT","THCA","UCEC")

data <- data.frame(
  Drinking = rep(c("Never", "Previous", 
                   "Current", "Dependent"), length(TCGA_CANCER)),
  Cancer = rep(TCGA_CANCER, each = 4),
  OR = NA,
  Pval = NA,
  Pval_FDR_cancer = NA_real_,   # 按癌种BH-FDR校正值
  stringsAsFactors = FALSE
)
ivw_need<-data.frame()
for (cancer_code in TCGA_CANCER) {
  cancer_info <- cancer_list[[cancer_code]]
  cancer_name <- cancer_info$name
  cancer_ids <- cancer_info$ids
  
  cancer_data <- res_ivw[res_ivw$id.outcome %in% cancer_ids, ]
  cancer_data <- merge(res_group, cancer_data, by = "id.exposure")
  
  drinking_groups <- c("Never", "Previous", 
                       "Current", "Dependent")
  cancer_data <- cancer_data[cancer_data$group %in% drinking_groups, ]
  cancer_data$group <- factor(cancer_data$group, levels = drinking_groups)
  #可视化
  cancer_data <- cancer_data[order(cancer_data$group), ]
  f <- function(x) {
    if(x < 0.01) {
      sprintf("%.2e", x)
    } else {
      round(x, 2)
    }
  }
  cancer_data$pval_num <- as.numeric(cancer_data$pval)
  cancer_data$pval <- sapply(cancer_data$pval, f)
  cancer_data$combined <- paste(round(cancer_data$b, 4), "(", cancer_data$pval, ")", sep = " ")
  prelim_meta <- metagen(
    TE = log(cancer_data$or),
    seTE = (log(cancer_data$or_uci95) - log(cancer_data$or_lci95)) / (2 * 1.96),
    studlab = paste(cancer_data$id.exposure, cancer_data$id.outcome, sep = "→"),
    data = cancer_data,
    sm = "OR",
    common = TRUE,
    random = TRUE
  )
  # 判断是否使用随机效应模型
  use_random <- prelim_meta$pval.Q < 0.1 | prelim_meta$I2 > 50
  
  
  meta_final <- metagen(
    TE = log(cancer_data$or),
    seTE = (log(cancer_data$or_uci95) - log(cancer_data$or_lci95)) / (2 * 1.96),
    studlab = paste(cancer_data$id.exposure, cancer_data$id.outcome, sep = "→"),
    data = cancer_data,
    sm = "OR",
    common = !use_random,
    random = use_random,
    subgroup = group,
    print.subgroup.name = FALSE,
    test.subgroup = FALSE
  )
  
  for (drinking_group in drinking_groups) {
    group_data <- cancer_data[cancer_data$group == drinking_group, ]
    if (nrow(group_data) == 0) {
      or <- 1.000
      pval <- NA
      pval_num <- NA    #无研究时FDR用NA
    }else{
      if (nrow(group_data) == 1) {
        or <- group_data$or[1]
        pval <- group_data$pval[1]
        pval_num <- group_data$pval_num[1]   
      } else {
        meta_final <- metagen(
          TE = log(group_data$or),
          seTE = (log(group_data$or_uci95) - log(group_data$or_lci95)) / (2 * 1.96),
          studlab = paste(group_data$id.exposure, group_data$id.outcome, sep = "→"),
          data = group_data,
          sm = "OR",
          common = !use_random,
          random = use_random
        )
        if (use_random) {
          or <- exp(meta_final$TE.random)
          pval <- meta_final$pval.random
        } else {
          or <- exp(meta_final$TE.common)
          pval <- meta_final$pval.common
        }
      }
      pval_num <- suppressWarnings(as.numeric(as.character(pval)))
      meta_summary <- data.frame(
        cancer_code = cancer_code,
        cancer_name = cancer_name,
        drinking_group = drinking_group,
        n_studies = nrow(group_data),
        OR = or,
        Pval = pval,
        I2 = meta_final$I2,
        Q_pvalue = meta_final$pval.Q,
        model = ifelse(use_random, "Random", "Fixed"),
        stringsAsFactors = FALSE
      )
      meta_results <- rbind(meta_results, meta_summary)
    }
    data_idx <- which(data$Cancer == cancer_code & data$Drinking == drinking_group)
    if (length(data_idx) > 0) {
      data$OR[data_idx] <- or
      data$Pval[data_idx] <- pval
      data$Pval_raw[data_idx] <- pval_num   
    }
  }
  
  

  # meta P值做 Benjamini-Hochberg FDR 校正，仅新增一列，不改变原OR/P计算。
  idx_cn <- which(data$Cancer == cancer_code)
  pv <- data$Pval_raw[idx_cn]   
  ok <- !is.na(pv)
  q <- rep(NA_real_, length(pv))
  if (sum(ok) > 0) {
    po <- pv[ok]
    ord <- order(po)
    n <- length(po)
    qo <- po[ord] * n / (1:n)
    qo <- pmin(rev(cummin(rev(qo))), 1)
    q[ok][ord] <- qo
  }
  data$Pval_FDR_cancer[idx_cn] <- q
  
  qg <- data$Pval_FDR_cancer[idx_cn]
  names(qg) <- data$Drinking[idx_cn]
  cancer_data$group_lab <- paste0(cancer_data$group,
                                  ifelse(!is.na(qg[cancer_data$group]) & qg[cancer_data$group] < 0.05, " *", ""))
  meta_final <- metagen(
    TE = log(cancer_data$or),
    seTE = (log(cancer_data$or_uci95) - log(cancer_data$or_lci95)) / (2 * 1.96),
    studlab = paste(cancer_data$id.exposure, cancer_data$id.outcome, sep = "→"),
    data = cancer_data,
    sm = "OR",
    common = !use_random,
    random = use_random,
    subgroup = group_lab,
    print.subgroup.name = FALSE,
    test.subgroup = FALSE
  )
  ivw_need<-rbind(ivw_need,cancer_data)
  pdf_file <- paste0("alcohol drinkers status/02. meta/", cancer_code, "_", gsub(" ", "_", cancer_name), "_meta_analysis.pdf")
  pdf(file = pdf_file, width = 14, height = 8)
  forest(meta_final, 
         col.diamond = "blue2",
         col.diamond.lines = "blue2",
         prediction = TRUE,
         print.tau2 = FALSE,
         leftcols = c("studlab", "nsnp", "combined"),
         leftlabs = c(paste("Study on alcohol drinker status and", cancer_code), "SNPs", "Beta (P value) (*: per-cancer FDR<0.05)"))
  
  dev.off()
}

write_xlsx(ivw_need, "alcohol drinkers status/ivw_results_FDR_drinking.xlsx")

