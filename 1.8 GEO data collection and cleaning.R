library(tidyverse)
library(limma) 
# GSE30219.Rdata was downloaded by GEOquery package
load("dat_input/GSE30219.Rdata")
platformMap <- data.table::fread("dat_input/platformMap.txt",data.table = F)

#1.7 Chip data analysis


exprSet <- exprs(get(ls(pattern = ".dat"))) 

ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), 
                          na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

## Start judging
if (LogC) { 
  ex[which(ex <= 0)] <- NaN
  ## Take log2
  exprSet <- log2(ex)
  print("log2 transform finished")
}else{
  print("log2 transform not needed")
}


#boxplot(exprSet,outline=FALSE, notch=T, las=2)
### By default, this function uses quntile to correct differences

exprSet=normalizeBetweenArrays(exprSet)
#boxplot(exprSet,outline=FALSE, notch=T, las=2)
## This step is important to convert the matrix into a data box
exprSet <- as.data.frame(exprSet)

#################################################################
## Probe gene name conversion
##platformMap There is a common platform R annotation package correspondence, this is what I put together。

index <- unique(get(ls(pattern = "pheno"))$platform_id)

 
if(any(grepl(index,platformMap$gpl))){


package <- paste0(platformMap$bioc_package[grep(index,platformMap$gpl)],".db")
## Install R package, can be installed directly, here use the judgment

if(!requireNamespace(package )){
  options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
  BiocManager::install(package ,update = F,ask = F)
} 

## Load R package
library(package,character.only =T)

## Obtaining the correspondence between the probe and the gene: This is a key step in the annotation of the probe
probe2symbol_df <- toTable(get(paste0(str_remove(package,".db"),"SYMBOL")))
dat_pre <- exprSet %>% 
  ## Row name to column name, because inner_join is only used for columns that become data boxes
  rownames_to_column("probe_id") %>% 
  ## Merge probe information
  inner_join(probe2symbol_df,by="probe_id") %>% 
  ## Remove superfluous information
  dplyr::select(-probe_id) %>%  
  ## rearrange
  dplyr::select(symbol,everything()) %>%  
  ## rowMeans average the number of trips. Represents the data passed in above)
  ## .[,-1] means that the first column in and out of the data is removed and the average of the rows is calculated
  mutate(rowMean =rowMeans(.[,-1])) %>% 
  ## Order the mean of the expressions from greatest to smallest
  arrange(desc(rowMean)) %>% 
  ## De-weight, symbol keeps the first one
  distinct(symbol,.keep_all = T) %>% 
  ## Reverse selection removes the rowMean column
  dplyr::select(-rowMean) %>% 
  ## Column name changes to row name
  column_to_rownames("symbol")

assign(paste0("GSE30219",".exp"),dat_pre)

dat.save <- ls(pattern = ".exp|pheno")

save(list = dat.save,
     file = paste0("dat_output/","GSE30219",".clean.Rdata"))

rm(list = ls(pattern = 'GSE'))

gc()

}

# View the result --------------------------------------------------------------------

## OS
load("dat_output/GSE30219.clean.Rdata")











