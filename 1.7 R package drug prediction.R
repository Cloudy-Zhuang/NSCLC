## package
rm(list = ls())  
library(tidyverse)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)

## Prepare data
load("dat_input/tcga_panmRNA_expr.Rdata")
load("dat_input/NSCLC-exp-clin-data.Rdata")
load("dat_output/通路基因生存分析及表达谱分析结果.Rdata")
ts <- tcga_panmRNA_expr[1:20]

nsclc_fpkm <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)
gene.important <- sort(setdiff(gene.important,"RPS6KA3"))

tcga_panmRNA_expr %>% 
  select(any_of(substr(names(nsclc_fpkm),1,15))) %>% 
  .[gene.important,] %>% 
  mutate_all(~.x^2-0.001)->
  tcga_clean



# Drug prediction --------------------------------------------------------------------


th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
dir='dat_input/'
# Representation matrix of training set
CTRP2_Expr = readRDS(file=file.path(dir,'CTRP2_Expr (TPM, not log transformed).rds'))
CTRP2_Res = readRDS(file = file.path(dir,"CTRP2_Res.rds"))
CTRP2_Res <- exp(CTRP2_Res) 


# A representation matrix that requires prediction ---------------------------------------------------------------
tcga_clean <- as.matrix(tcga_clean)


calcPhenotype(trainingExprData = CTRP2_Expr ,
              trainingPtype = CTRP2_Res,
              testExprData =  tcga_clean,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

