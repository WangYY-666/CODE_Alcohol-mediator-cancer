
rm(list=ls())
library(dplyr)

load(file = "prime/2/diffgenes.Rdata")
load(file = "ESCA_tpm_clin (176 samples).Rdata")
#load(file = "clinic_data(types).Rdata")
gene<-gene_name$Gene

#tpm
#load(file = "diffgenes_CA.Rdata")
#gene<-gene_name$allg
rt_cli<-rt_tpm[rt_tpm$sample %in% clin_CA$sample,]
rt_cli1<-rt_cli[,gene]
rt_cli2<-cbind(rt_cli[,c(11,12)],rt_cli1)
housing.df<-as.data.frame(lapply(rt_cli2,as.numeric))
rownames(housing.df)<-rt_cli$id

#count
load(file = "count_ESCA(clear).Rdata")
rt_ESCA<-expr_ESCA
rt1<-rt_ESCA[gene,]
rt1<-log2(rt1+1)
rt2<-as.data.frame(t(rt1))
rt2<-cbind(sample=rownames(rt2),rt2)
rt3<-merge(rt_cli[,c(1:13)],rt2,by="sample")
rt4<-rt3[,c(2,11,12,14:ncol(rt3))]
housing.df<-column_to_rownames(rt4,var = "id")


with_dash <- grep("-", gene, value = TRUE)
#add<-cbind(rt_cli[,c(11,12)],rt_cli1[,with_dash])
add<-cbind(rt4[,c(2,3)],rt4[,with_dash])
with_dash1<-sub("-.*", "", with_dash)
colnames(add)[3:6]<-with_dash1<-sub("-.*", "", colnames(add)[3:6])
univ_formulas <- sapply(with_dash1,
                        function(x) as.formula(paste('Surv(OS.time,OS)~', x)))
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = add)})



without_dash <- grep("-", gene, value = FALSE, invert = TRUE)
without_dash <- gene[without_dash]  

univ_formulas <- sapply(without_dash,
                        function(x) as.formula(paste('Surv(OS.time,OS)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = housing.df)})

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)

                         p.value<-signif(x$wald["pvalue"], digits=2)
                      
                         HR <-signif(x$coef[2], digits=2);
                     
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })

res <- t(as.data.frame(univ_results, check.names = FALSE))
res<-as.data.frame(res)

#res_p<-res[res$p.value<0.05,]

rownames(res)<-with_dash
res_p1<-res[res$p.value<0.05,]

res_p2<-res[res$p.value<0.05,]

res_p<-rbind(res_p1,res_p2)

save(res_p,file = "LASSO/res_p.Rdata")

housing.df1<-housing.df[,rownames(res_p)]
housing.df2<-housing.df[,c(1,2)]
housing.df<-cbind(housing.df2,housing.df1)
housing.df$OS.time<-housing.df$OS.time/365

colnames(housing.df)[3]<-"MIR7_3HG"

set.seed(123) ## to get the same sequence of numbers
# randomly sample 60% of the row IDs for training; the remaining 40% serve as validation
train.rows <- sample(rownames(housing.df), dim(housing.df)[1]*0.6)

# collect all the columns with training row ID into training set:
housing.train <- housing.df[train.rows, ]
save(housing.train,file = "LASSO/traindata.Rdata")
#housing.train <- housing.df

valid.rows <- setdiff(rownames(housing.df), train.rows)
housing.valid <- housing.df[valid.rows, ]
save(housing.valid,file = "LASSO/housing_valid.Rdata")
#write.table(housing.valid,"housing.valid.txt",sep = "\t",quote = F)
glmnet_input=housing.train
glmnet_input<-glmnet_input[order(glmnet_input[,1]),]

x <- data.matrix(glmnet_input[,3:length(colnames(glmnet_input))])

#最佳alpha
#library(caret)
#model <- train(OS ~., data = glmnet_input,method = 'glmnet', trControl = trainControl('cv', number = 10),  tuneLength = 10)
#model$bestTune  #最佳 alpha 和 λ 值
#coef(model$finalModel, model$bestTune$lambda)  #各预测变量的回归系数

library(survival)
y <- data.matrix(Surv(glmnet_input[,1],glmnet_input[,2]))
## lasso
library(glmnet)
cv.fit <- cv.glmnet(x, y,family="cox") 
#cv.fit<-cv.glmnet(x,y,family="cox", alpha=0.5)
cv.fit$lambda.min
cv.fit$lambda.1se
save(cv.fit,file = "LASSO/cv.fit.Rdata")

#x11()
plot(cv.fit)

if(T){
  library(ggplot2)
  library(broom)
  tidied_cv = tidy(cv.fit)
  ggplot(tidied_cv, aes(lambda, estimate)) +
    geom_point(color = "red",size=2) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), alpha = .3) +
    scale_x_log10()+
    geom_vline(xintercept = cv.fit$lambda.min, lty=2,color = "red",size=1) +
    geom_vline(xintercept = cv.fit$lambda.1se, lty=2) +
    ylab("Partial Likelihood Deviance")+
    annotate("text",x=cv.fit$lambda.min/5,size=5,
             y =max(tidied_cv$estimate),
             label =paste0("lambda.min = ",
                           round(cv.fit$lambda.min,3),
                           "\n",
                           "number = ",tidied_cv$nzero[tidied_cv$lambda==cv.fit$lambda.min]))+
    ggtitle("")+
    theme_bw()
}

coef1<-coef(cv.fit,s = "lambda.min")
coef2<-coef(cv.fit,s = "lambda.1se")
coef1
coef2
#fit <- glmnet(x, y, family = "cox", alpha=1)
fit <- glmnet(x, y, family = "cox")
fit_predict <- fit
save(fit_predict,file = "LASSO/fit_predict.Rdata")

library(broom)
library(patchwork) 
library(tidyverse)
tidy_df <- broom::tidy(fit)            
tidy_cvdf <- broom::tidy(cv.fit)            
tmp <- filter(tidy_df, term != '(Intercept)')            
ggplot(tmp, aes(log(lambda),estimate,color=term,group=term))+     geom_line(size=1)+            
  labs(x="Log Lambda",y="Coefficients")+            
  coord_cartesian(xlim = c(-9,-1), expand = T) +            
  guides(color = guide_legend(title = ''),direction = 'cox') +            
  #scale_x_continuous(limits = c(-9,-1)) +            
  theme_bw()

plot(fit, xvar = "lambda", label = TRUE)
if(T){
  (Coefficients <- coef(fit, s = cv.fit$lambda.min))
  (Active.Index <- which(as.matrix(Coefficients) != 0))
  (Active.Coefficients  <- Coefficients[Active.Index])
  
  (gene_found <- row.names(Coefficients)[Active.Index])
  save(gene_found,file = "LASSO/gene_found.Rdata")
}

if(T){
  tidy_covariates = tidy(cv.fit$glmnet.fit)
  index <- tidy_covariates$term %in% gene_found
  library(dplyr)
  data_point <- tidy_covariates %>% 
    filter(index) %>% 
    arrange(lambda) %>% 
    distinct(term,.keep_all = T)
  library(ggrepel)
  library(ggplot2)
  ggplot(tidy_covariates)+
    geom_line(data = subset(tidy_covariates,index), aes(x=lambda, y=estimate, color=as.factor(term)),size=1) +
    geom_line(data = subset(tidy_covariates,!index), aes(x=lambda, y=estimate, color=as.factor(term)),size=0.8,alpha=0.1) +
    guides(color=FALSE) +
    geom_vline(xintercept = cv.fit$lambda.min, lty=2,color='red') +
    scale_x_log10()+
    ylab("Coefficients")+
    geom_point(data =data_point,aes(lambda,estimate,color=as.factor(term)),size=3)+
    geom_text_repel(data =data_point,aes(lambda,estimate,label=term),size=3)+
    theme_classic()
}
if(T){
  coxdata <- cbind(glmnet_input[,c(1,2)],glmnet_input[,gene_found])
  res <- data.frame()
  genes = gene_found
  for (i in 1:length(genes)) {
    print(i)
    surv =as.formula(paste('Surv(OS.time,OS)~', genes[i]))
    cox = coxph(surv, data = coxdata)
    cox = summary(cox)
    p.value=signif(cox$wald["pvalue"], digits=2)
    HR =signif(cox$coef[2], digits=2);#exp(beta)
    HR.confint.lower = signif(cox$conf.int[,"lower .95"], 2)
    HR.confint.upper = signif(cox$conf.int[,"upper .95"],2)
    CI <- paste0("(", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    res[i,1] = genes[i]
    res[i,2] = HR
    res[i,3] = CI
    res[i,4] = p.value
  }
  names(res) <- c("ID","HR","95% CI","p.value")
  save(res,file = "LASSO/res_HR.Rdata")
}


coxresult <- cbind(res[,],Active.Coefficients)

colnames(coxresult) <- c("Symbol", 
                         "HR", "95% CI", "p.value", "Cox Coefficient")
write.table(coxresult,file = "LASSO/coefficience.txt",sep = "\t",row.names = F,quote = F)
coxresult$HR <- as.numeric(as.character(coxresult$HR))
coxresult <-  dplyr::arrange(coxresult,desc(HR))


x <- data.matrix(glmnet_input[,3:length(colnames(glmnet_input))])

rt <- glmnet_input  
rt <- rt[,1:2]
order =rownames(rt)
library(glmnet)
rt$riskScore <- predict(fit_predict, x, s=cv.fit$lambda.min, type="link")
# best AUC cutoff
rt <- as.data.frame(apply(rt,2,as.numeric))
rownames(rt) <- rownames(glmnet_input)

# rt1<-cbind(id=rownames(rt),rt)
# rt2<-cbind(id=rownames(glmnet_input),glmnet_input[,c(3:6)])
#rt3<-merge(rt1,rt2,by="id")
#save(rt3,file = "cox_riskscore.Rdata")


source("survivalROC_NEW.R")
cutoffROC <- function(t,rt){
  rt=rt
  ROC_rt <- survivalROC_NEW(Stime=rt[,1], status=rt[,2], marker = rt[,3], 
                            predict.time =t, method="KM")
  
  youdenindex=max(ROC_rt$sensitivity+ROC_rt$specificity-1)
  
  index = which(ROC_rt$sensitivity+ROC_rt$specificity-1==youdenindex)+1
  
  result=c(ROC_rt$cut.values[index],ROC_rt$sensitivity[index-1],ROC_rt$specificity[index-1])
}
#t<-round(sum(housing.df$OS.time)/nrow(housing.df))
t<-sum(housing.train$OS.time)/nrow(housing.train)
#t <- 2 #选定预测时间的最佳截点
(res.cut <- cutoffROC(t,rt))
cutoff <- res.cut[1]



library(dplyr)
if(T){
  rt$risk=as.vector(ifelse(rt$riskScore>cutoff,"high","low"))
  rt$order <- order
  rt$riskScore = as.numeric(rt$riskScore)
  rt = rt %>% arrange(riskScore)
  rt$num =seq(1,length(rownames(rt)))
}
rt1<-rt
####
#save(rt1,file = "18need.Rdata") 

if(T){
  cutnum <- sum(rt$risk=="low")
  library(ggplot2)
  (plot.point=ggplot(rt,aes(x=num,y=riskScore))+
      geom_point(aes(col=risk),size=3)+
      geom_segment(aes(x = cutnum, y = min(riskScore), xend = cutnum, yend = cutoff),linetype="dashed")+
      geom_segment(aes(x = 0, y = cutoff, xend = cutnum,yend = cutoff),linetype="dashed")+
      geom_text(aes(x=cutnum/2,y=cutoff+0.25,label=paste0("Cutoff: ",round(cutoff,2))),col ="black",size = 8,alpha=0.8)+
      theme_classic()+
      theme(axis.title.x=element_blank(),#legend.text = element_text(size = 15),
            axis.text.x=element_text(size = 20),
            #axis.ticks.x=element_blank(),
            axis.title.y =  element_text(size = 20),
            axis.text.y = element_text(size = 15))+
      scale_color_manual(values=c("#FC4E07","#00AFBB"))+
      scale_x_continuous(limits = c(0,NA),expand = c(0,0))
  )
}


if(T){
  event = factor(rt$OS)
  (plot.sur=ggplot(rt,aes(x=num,y=OS.time))+
      geom_point(aes(col=event),size=3)+
      geom_vline(aes(xintercept = cutnum),linetype="dashed")+
      theme_classic()+
      theme(axis.title.x=element_text(size = 15),
            axis.text.x=element_text(size = 20),
            axis.ticks.x=element_blank(),
            axis.title.y =  element_text(size = 20),
            axis.text.y = element_text(size = 15))+      
      scale_color_manual(values=c("#00AFBB","#FC4E07"))+
      scale_x_continuous(limits = c(0,NA),expand = c(0,0)))
}


if(T){
  tmp=t(glmnet_input[rt$order,rev(coxresult$Symbol)])
  library(data.table)
  library(plyr)
  library(scales)
  tmp.m <- reshape2::melt(tmp)
  colnames(tmp.m)=c("Name", "variable", "value")
  tmp.m <- ddply(tmp.m, .(variable), transform,rescale = rescale(value))
  (plot.h <- ggplot(tmp.m, aes(variable, Name)) + 
      geom_tile(aes(fill = rescale), colour = "white") + 
      scale_fill_gradient(low = "#00AFBB",high = "#FC4E07")+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y = element_text(size = 15))
  )
}

library(cowplot)
plot_grid(plot.point, plot.sur, plot.h,
          align = 'v',ncol = 1,axis="b")

##########################################################################################



if(T){
  library("survival")
  library("survminer")
  my.surv <- Surv(rt$OS.time, rt$OS)
  group <- rt$risk
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  
  
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  x = summary(coxph(Surv(OS.time, OS)~riskScore, data = rt))
  HR = signif(x$coef[2], digits=2)
  up95 = signif(x$conf.int[,"upper .95"],2)
  low95 = signif(x$conf.int[,"lower .95"], 2)
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  rt <- rt[order(rt[,"riskScore"],decreasing = T),]
  plot.km<-ggsurvplot(fit, data = survival_dat ,
                      
                      conf.int = F, 
                      
                      censor = F,
                      palette = c("#D95F02","#1B9E77"), 
                      
                      font.legend = 11,
                      
                      legend.labs=c(paste0(">",round(cutoff,2),"(",sum(rt$risk=="high"),")"),
                                    paste0("<",round(cutoff,2),"(",sum(rt$risk=="low"),")")),
                      ylab=c('Survival probability'),
                      
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
}


##ROC
if(T){
  library(tidyr)
  library(purrr)
  library(survivalROC)
  ## Define a helper function to evaluate at various t
  survivalROC_helper <- function(t) {
    survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$riskScore, 
                predict.time =t, method="KM")
  }
  ## Evaluate 3,5,10
  survivalROC_data <- data_frame(t = c(1,2,3)) %>%
    mutate(survivalROC = map(t, survivalROC_helper),
           ## Extract scalar AUC
           auc = purrr::map_dbl(survivalROC, magrittr::extract2, "AUC"),
           ## Put cut off dependent values in a data_frame
           df_survivalROC = map(survivalROC, function(obj) {
             dplyr::as_data_frame(obj[c("cut.values","TP","FP")])
           })) %>%
    dplyr::select(-survivalROC) %>%
    unnest() %>%
    arrange(t, FP, TP)
  ## Plot
  survivalROC_data1 <- survivalROC_data %>% 
    mutate(auc =sprintf("%.3f",auc))%>% 
    unite(year, t,auc,sep = " years AUC: ")
  
  AUC =factor(survivalROC_data1$year)
  plot.roc<-survivalROC_data1 %>%
    ggplot(mapping = aes(x = FP, y = TP)) +
    geom_path(aes(color= AUC),size=2)+
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    theme_classic() +
    theme(legend.position = c(0.7,0.2))+theme(axis.title.x=element_blank(),
                                              axis.text.x=element_blank(),
                                              axis.ticks.x=element_blank(),
                                              axis.title.y =  element_text(size = 15),
                                              axis.text.y = element_text(size = 15),
                                              legend.text = element_text(size = 12))
}


plot_grid(plot.km$plot, plot.roc,
          align = 'v',ncol = 1,axis="b")











#***validation
glmnet_input=housing.valid
x <- data.matrix(glmnet_input[,3:length(colnames(glmnet_input))])

rt <- glmnet_input  
rt <- rt[,1:2]
order =rownames(rt)

library(glmnet)
rt$riskScore <- predict(fit_predict, x, s=cv.fit$lambda.min, type="link")

rt <- as.data.frame(apply(rt,2,as.numeric))
rownames(rt) <- rownames(glmnet_input)

if(T){
  rt$risk=as.vector(ifelse(rt$riskScore>cutoff,"high","low"))
  rt$order <- order
  rt$riskScore = as.numeric(rt$riskScore)
  rt = rt %>% arrange(riskScore)
  rt$num =seq(1,length(rownames(rt)))
}
rt2<-rt
####

source("survivalROC_NEW.R")



# 准备数据
library(dplyr)
if(T){
  rt$risk=as.vector(ifelse(rt$riskScore>cutoff,"high","low"))
  rt$order <- order
  rt$riskScore = as.numeric(rt$riskScore)
  rt = rt %>% arrange(riskScore)
  rt$num =seq(1,length(rownames(rt)))
}


if(T){
  cutnum <- sum(rt$risk=="low")
  library(ggplot2)
  (plot.point=ggplot(rt,aes(x=num,y=riskScore))+
      geom_point(aes(col=risk),size=3)+
      geom_segment(aes(x = cutnum, y = min(riskScore), xend = cutnum, yend = cutoff),linetype="dashed")+
      geom_segment(aes(x = 0, y = cutoff, xend = cutnum,yend = cutoff),linetype="dashed")+
      geom_text(aes(x=cutnum/2,y=cutoff+0.25,label=paste0("Cutoff: ",round(cutoff,2))),col ="black",size =8,alpha=0.8)+
      theme_classic()+
      theme(axis.title.x=element_blank(),#legend.text = element_text(size = 15),
            axis.text.x=element_text(size = 20),
            #axis.ticks.x=element_blank(),
            axis.title.y =  element_text(size = 20),
            axis.text.y = element_text(size = 15))+
      scale_color_manual(values=c("#FC4E07","#00AFBB"))+
      scale_x_continuous(limits = c(0,NA),expand = c(0,0)))
}


if(T){
  event = factor(rt$OS)
  (plot.sur=ggplot(rt,aes(x=num,y=OS.time))+
      geom_point(aes(col=event),size=3)+
      geom_vline(aes(xintercept = cutnum),linetype="dashed")+
      theme_classic()+
      theme(axis.title.x=element_text(size = 15),
            axis.text.x=element_text(size = 20),
            axis.ticks.x=element_blank(),
            axis.title.y =  element_text(size = 20),
            axis.text.y = element_text(size = 15))+
      scale_color_manual(values=c("#00AFBB","#FC4E07"))+
      scale_x_continuous(limits = c(0,NA),expand = c(0,0)))
}


if(T){
  tmp=t(glmnet_input[rt$order,rev(coxresult$Symbol)])
  library(data.table)
  library(plyr)
  library(scales)
  tmp.m <- reshape2::melt(tmp)
  colnames(tmp.m)=c("Name", "variable", "value")
  tmp.m <- ddply(tmp.m, .(variable), transform,rescale = rescale(value))
  (plot.h <- ggplot(tmp.m, aes(variable, Name)) + 
      geom_tile(aes(fill = rescale), colour = "white") + 
      scale_fill_gradient(low = "#00AFBB",high = "#FC4E07")+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y = element_text(size = 15)))
}

library(cowplot)
plot_grid(plot.point, plot.sur, plot.h,
          
          align = 'v',ncol = 1,axis="b")
##########################################################################################



if(T){
  library("survival")
  library("survminer")
  my.surv <- Surv(rt$OS.time, rt$OS)
  group <- rt$risk
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  
  
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  x = summary(coxph(Surv(OS.time, OS)~riskScore, data = rt))
  HR = signif(x$coef[2], digits=2)
  up95 = signif(x$conf.int[,"upper .95"],2)
  low95 = signif(x$conf.int[,"lower .95"], 2)
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  rt <- rt[order(rt[,"riskScore"],decreasing = T),]
  plot.km<-ggsurvplot(fit, data = survival_dat ,
                      
                      conf.int = F, 
                      
                      censor = F,
                      palette = c("#D95F02","#1B9E77"), 
                      
                      font.legend = 11,
                      
                      legend.labs=c(paste0(">",round(cutoff,2),"(",sum(rt$risk=="high"),")"),
                                    paste0("<",round(cutoff,2),"(",sum(rt$risk=="low"),")")),
                      ylab=c('Survival probability'),
                      
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
}

##ROC
if(T){
  library(tidyr)
  library(purrr)
  library(survivalROC)
  ## Define a helper function to evaluate at various t
  survivalROC_helper <- function(t) {
    survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$riskScore, 
                predict.time =t, method="KM")
  }
  ## Evaluate 3,5,10
  survivalROC_data <- data_frame(t = c(1,2,3)) %>%
    mutate(survivalROC = map(t, survivalROC_helper),
           ## Extract scalar AUC
           auc = purrr::map_dbl(survivalROC, magrittr::extract2, "AUC"),
           ## Put cut off dependent values in a data_frame
           df_survivalROC = map(survivalROC, function(obj) {
             dplyr::as_data_frame(obj[c("cut.values","TP","FP")])
           })) %>%
    dplyr::select(-survivalROC) %>%
    unnest() %>%
    arrange(t, FP, TP)
  ## Plot
  survivalROC_data1 <- survivalROC_data %>% 
    mutate(auc =sprintf("%.3f",auc))%>% 
    unite(year, t,auc,sep = " years AUC: ")
  
  AUC =factor(survivalROC_data1$year)
  plot.roc<-survivalROC_data1 %>%
    ggplot(mapping = aes(x = FP, y = TP)) +
    geom_path(aes(color= AUC),size=2)+
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    theme_classic() +
    theme(legend.position = c(0.7,0.2))+theme(axis.title.x=element_blank(),
                                              axis.text.x=element_blank(),
                                              axis.ticks.x=element_blank(),
                                              axis.title.y =  element_text(size = 15),
                                              axis.text.y = element_text(size = 15),
                                              legend.text = element_text(size = 12))
}

plot_grid(plot.km$plot, plot.roc,
          align = 'v',ncol = 1,axis="b")

