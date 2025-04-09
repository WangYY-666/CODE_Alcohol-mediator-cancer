
rm(list=ls())

library(IOBR)
library(magrittr)
library(ggh4x)
library(tidyverse)
library(readr)
library(patchwork)
library(stringr)
########################################################

cols <- IOBR::palettes(category = "random", palette = 4,
                       show_col = F,
                       show_message = F)
cib_exist_all= exist_cib%>%
  subset(.,CancerType=='ESCA')%>%#取ESCA
  dplyr::rename(sample=SampleID)%>%
  mutate(sample=gsub('\\.','-',sample))%>%
  dplyr::select(-CancerType)%>% #去掉这一列
  dplyr::select(-c((ncol(.) - 2):ncol(.)))%>%#去掉后三列
  na.omit()%>%
  distinct(sample,.keep_all = T)
#cib_exist_ESCA<-cibersort[,-c(24:26)]
##
load(file = "Clinic data(ESCA). Rdata")
load(file = "id_ESCA_all.Rdata")

#EAC
cib_exist_CA<-merge(id_CA[,c(1,3)],cib_exist_all,by="sample")
cib_exist_ESCA<-cib_exist_CA
#SCC
cib_exist_SCC<-merge(id_SCC[,c(1,3)],cib_exist_all,by="sample")
cib_exist_ESCA<-cib_exist_SCC
meltdf_cib_exist = cib_exist_ESCA%>%
  reshape2::melt()%>%
  mutate(sample=factor(sample,levels=cib_exist_ESCA$sample))%>%
  mutate(Alcohol_history=factor(alcohol_history,levels=c('No','Yes')))

plot_tme = function(meltdf){
  
  
  p_stack_df = ggplot(meltdf,aes(x = sample, y = value, fill = variable)) + 
    geom_bar(stat = "identity")+ theme_test() + 
    scale_fill_manual(values = cols) + #scale_x_discrete(limits = rev(levels(input))) + 
    #ggtitle(paste0(title)) + 
    theme(legend.position = 'bottom')+
    theme(plot.title = element_text(size = rel(2),hjust = 0.5), 
          axis.text.x = element_blank())+
    theme(axis.ticks =element_blank() )+
    scale_y_continuous(expand = c(0,0))+
    facet_nested(.~Alcohol_history,drop=T,scale="free",space="free",switch="x",
                 strip =strip_nested(background_x = elem_list_rect(fill =c('steelblue','firebrick')),by_layer_x = F))
  
  # p_stack_df
  
  p_box_df = ggplot(meltdf,aes(x = variable, y = value, 
                               fill = variable)) +
    geom_boxplot(width=0.4,lwd=0.2,color='black',outlier.shape=NA,
                 position = position_dodge(width = 0.8))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5))+
    theme(legend.position = 'none')+
    scale_fill_manual(values = cols)
  # p_box_df
  
  p_box_comp =  ggplot(meltdf,aes(x = variable, y = value, 
                                  fill = Alcohol_history,color = Alcohol_history)) +
    
    geom_boxplot(width=0.5,lwd=0.2,color='black',outlier.shape=NA,
                 position = position_dodge(width = 0.8))+
    theme_bw()+
    theme(legend.position = 'top')+
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5))+
    stat_compare_means(aes(group=Alcohol_history,label = after_stat(p.signif)),hide.ns = FALSE,
                       method = "t.test")+
    scale_fill_manual(values = c('steelblue','firebrick'))
  
  # p_box_comp
  
  
  p_heat = ggplot(meltdf,aes(x = sample, y =variable ,fill=value))+
    geom_tile()+
    scale_fill_gradientn(colors = c('white','steelblue','firebrick'))+
    
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(fill = 'NA',color = 'black'))+
    labs(fill="Z score", x = NULL, y = NULL)
  # p_heat
  plot_comb = (p_stack_df)/(p_box_df|p_box_comp)/p_heat
  return(plot_comb)
}
plot_exist_cib = plot_tme(meltdf_cib_exist)
pdf('./plot_cibersort_ESCA.pdf',width = 14,height = 14)
print(plot_exist_cib)
dev.off()


load(file = "cibersort_tpm.Rdata")
cibersort<-cibersort_tpm
colnames(cibersort)[1]<-"sample"
cib_exist_ESCA1<-merge(id_CA[,c(1,3)],cibersort[,1:23],by="sample")
cib_exist_ESCA1<-merge(id_SCC[,c(1,3)],cibersort[,1:23],by="sample")
meltdf_cib_exist1 = cib_exist_ESCA1%>%
  reshape2::melt()%>%
  mutate(sample=factor(sample,levels=cib_exist_ESCA1$sample))%>%
  mutate(Alcohol_history=factor(alcohol_history,levels=c('No','Yes')))
meltdf_cib_exist1$variable<-str_remove(meltdf_cib_exist1$variable, "_CIBERSORT$")
order1<-c("B_cells_memory","B_cells_naive" ,"Plasma_cells"   ,
          "T_cells_CD4_memory_activated", "T_cells_CD4_memory_resting" , "T_cells_CD4_naive"  ,         
          "T_cells_CD8" , "T_cells_follicular_helper" , "T_cells_gamma_delta",
          "T_cells_regulatory_(Tregs)","NK_cells_activated", "NK_cells_resting",
          "Macrophages_M0","Macrophages_M1" ,"Macrophages_M2",
          "Dendritic_cells_activated" , "Dendritic_cells_resting", "Mast_cells_activated",   
          "Mast_cells_resting" , "Monocytes" ,"Eosinophils" ,"Neutrophils")
meltdf_cib_exist1<-meltdf_cib_exist1[order(match(meltdf_cib_exist1$variable, order1)), ]
#meltdf_cib_exist1$variable<-gsub("_", ".", meltdf_cib_exist1$variable)
plot_exist_cib1 = plot_tme(meltdf_cib_exist1)

pdf('./plot_cibersort_ESCA(tpm).pdf',width = 14,height = 14)
print(plot_exist_cib1)
dev.off()

colnames(cibersort1)[1]<-"sample"
cib_exist_ESCA2<-merge(clin_ESCA[,c(1,13)],cibersort1[,1:23],by="sample")

meltdf_cib_exist2 = cib_exist_ESCA2%>%
  reshape2::melt()%>%
  mutate(sample=factor(sample,levels=cib_exist_ESCA2$sample))%>%
  mutate(Alcohol_history=factor(alcohol_history,levels=c('No','Yes')))
meltdf_cib_exist2$variable<-str_remove(meltdf_cib_exist2$variable, "_CIBERSORT$")
plot_exist_cib2 = plot_tme(meltdf_cib_exist2)

