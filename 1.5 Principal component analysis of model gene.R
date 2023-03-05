setwd("C:/Users/Cloudy/Desktop/突击文章/dat_out")
rm(list = ls())
# Load package
library(tidyverse)
library(factoextra)
library(FactoMineR)

# Read data
load("dat_input/NSCLC-exp-clin-data.Rdata")
gene.important <- c("ADCY2", "ARHGAP11A", "CAMK4", "CCNA2", "CCNB2", "CDCA5", "CENPH", "CFL1", "DSCC1", "FASLG", "MAD2L1", "MAPK10", "MKNK1", "PIK3C3", "PIK3R1", "PIK3R3", "PLK4", "PRKCB", "RAC1", "RAF1", "TPX2")

NSCLC_fpkm <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)

NSCLC_fpkm %>% 
  t() %>% 
  as.data.frame() %>% 
  select(gene.important) %>% 
  mutate(group=ifelse(substr(rownames(.),14,15)=="11","normal","cancer"))-> 
  pca_pre




tcga.pca <- prcomp(  pca_pre[,-22],scale = TRUE)


fviz_pca_ind(tcga.pca,
             label = "none", # hide individual labels
             habillage =   pca_pre$group, # color by groups
             palette = c("#00AFBB", "#E7B800","#FC4E07"),
             addEllipses = TRUE # Concentration ellipses
             )


# GEO Principal component analysis ---------------------------------------------------------------


load("GEO纳入研究的表达矩阵.Rdata")
load("NSCLC-COVID-drug单因素cox分析结果.Rdata")
## GSE14507

GSE147507_include_exp[cox_res$gene,] %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate_all(~log2(.x+1)) %>% 
  filter(!str_detect(rownames(.),"Lung")) %>% 
  mutate(group=ifelse(str_detect(rownames(.),"SARS-CoV|COVID19"),
                      "COVID-19","normal")) %>% 
  mutate(group=factor(group,levels = c("normal","COVID-19")))->
  GSE147507_dat




GSE147507_pca <- prcomp(GSE147507_dat[,-12],scale = F)


fviz_pca_ind(GSE147507_pca  ,
             label = "none", # hide individual labels
             habillage = GSE147507_dat$group, # color by groups
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE # Concentration ellipses
)

## It didn't work. I gave up



