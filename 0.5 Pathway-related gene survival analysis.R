rm(list = ls())
library(tidyverse)
library(survival)
library(survminer)

# bath uni-cox ------------------------------------------------------------
## read data

load("dat_input/NSCLC-exp-clin-data.Rdata")
load("dat_output/The PI3K AKT Mtor pathway genes for study. Rdata")


## Interest gene for previous literature
interest_gene <- sort(unique(c(mtor_hubgene,motor_gene)))
interest_gene <- sort(replace(interest_gene,interest_gene=="GRK2","ADRBK1"))

# Pathway-related genes in the literature
gene_publish <- c("PIK3R1", "PIK3CA", 
                  "PIK3C2A", "PIK3C3", 
                  "AKT1", "PIK3R2", "PIK3CB", 
                  "PIK3C2B", "PIK3R4", "AKT2",
                  "PIK3R3", 
                  "PIK3CD", "PIK3C2G", "PTEN", 
                  "AKT3", "TSC2", "MTOR","RHEB",
                  "PDK1","INPP4B","PHLPP1","STK11",
                  "TSC1")

## [Additional GENE]from liternatures FRAP1 SPP1、EFNA2、ITGA6、ITGB8、ITGB4、IL2RB、TNC、KLF5、HRH3
gene_new <- c("MTOR" , "SPP1" , "EFNA2", "ITGA6" ,"ITGB8",
              "ITGB4" ,"IL2RB" ,"TNC"  , "HRH3" )
gene_publish <- c(gene_publish,gene_new)
gene_publish <- sort(gene_publish)

gene.all <- sort(unique(c(interest_gene,gene_publish)))

mygene <- gene.all 
# create the NSCLC log2(fpkm+1) expression matrix
NSCLC_fpkmanno <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)
NSCLC_clin_tab <- rbind(LUAD_clin_tab,LUSC_clin_tab)

## import function for survival analysis
options(digits = 3)	
if (T) {
  #simple_sur
  
  simple_sur <- function(gene,dat){
    print(gene)
    #data=sur_dat_fin
    data_for_surv=select(dat,all_of(gene),1,2)
    # Calculate the optimal classification point
    x <-try(surv_cutpoint(data_for_surv,
                          time = "surv_time", 
                          event = "surv_status", 
                          variables = gene,
                          minprop = 0.25), silent=TRUE)
    if ('try-error' %in% class(x)){
      
      surv_res = matrix(c(gene,"无法分组",rep(NA,8)), nrow = 1)
      colnames(surv_res) = c("gene","cutoff_value","Num_Low", "Num_High", "Log-rank p", "coef", "HR", "Cox p", "Lower95%", "Upper95%")
      return(as_tibble(surv_res))
      
    }else{
      
      surv.cut = surv_cutpoint(data_for_surv,
                               time = "surv_time", 
                               event = "surv_status", 
                               variables = gene,
                               minprop = 0.25)}
    
    
    # Optimal classification value
    opticut = surv.cut$cutpoint$cutpoint 
    data_for_surv$group = surv_categorize(surv.cut)[,gene]
    
    data_for_surv$group=factor( data_for_surv$group,levels = c("low","high"))
    
    # Another situation that cannot be grouped
    if (any(is.na(data_for_surv$group))){
      
      surv_res = matrix(c(gene,"无法分组",rep(NA,8)), nrow = 1)
      colnames(surv_res) = c("gene","cutoff_value","Num_Low", "Num_High", "Log-rank p", "coef", "HR", "Cox p", "Lower95%", "Upper95%")
      return(as_tibble(surv_res))}
    
    
    # survival analysis
    # survfit object
    surv.obj = as.formula(paste0("Surv(surv_time, surv_status)~","group"))
    # K-M survival curves
    fit = surv_fit(surv.obj, data = data_for_surv) # 用survival::survfit报错: https://github.com/kassambara/survminer/issues/403
    # log-rank p http://rpkgs.datanovia.com/survminer/reference/surv_pvalue.html
    pval = surv_pvalue(fit)$pval
    # Cox proportional hazards regression model
    cox = coxph(surv.obj, data = data_for_surv)
    cox_summary = summary(cox)
    cox_res = c(cox_summary$coefficients[,c(1:2,5)],
                cox_summary$conf.int[,c(3,4)])
    
    high_pos <- paste(which(data_for_surv$group%in%"high"),collapse = "_")
    low_pos <-  paste(which(data_for_surv$group%in%"low"),collapse = "_")
    
   # Extract result output
    surv_res = matrix(c(gene,opticut,as.numeric(table(data_for_surv$group)), pval, cox_res,high_pos,low_pos), nrow = 1)
    colnames(surv_res) = c("gene","cutoff_value","Num_Low", "Num_High", "Log-rank p", "coef", "HR", "Cox p", "Lower95%", "Upper95%","high_pos",
                           "low_pos")
    # Extract result output
    as_tibble(surv_res)
 
  }
  
}


## create survival data
datClean <- function(fpkmanno,clin_tab,mygene){

mygene <- mygene[which(mygene%in%rownames(fpkmanno))] 


fpkmanno %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::mutate(id=substr(rownames(.),1,12)) %>% 
  dplyr::mutate(ind=rowSums(across(where(is.numeric)))) %>% 
  dplyr::arrange(desc(ind)) %>% 
  dplyr::distinct(id,.keep_all = T) %>% 
  dplyr::select(id,any_of(mygene)) %>% 
  dplyr::inner_join(clin_tab %>% 
               rownames_to_column("id")) %>% 
  column_to_rownames("id") %>% 
  dplyr::select(any_of(sort(mygene)),everything()) %>% 
  dplyr::rename("surv_time"=OS.time,"surv_status"=OS) %>% 
  dplyr::select(surv_time,surv_status,everything())->
    dat_clean 

  
  map_df(names(dat_clean)[3:which(names(dat_clean)==sort(mygene)[length(mygene)])],
         simple_sur,
         dat_clean) %>% 
    filter(`Log-rank p`<0.05&`Cox p`<0.05) %>% 
    mutate(across(-c(gene,high_pos,low_pos),as.numeric)) %>% 
    mutate(across(where(is.numeric),round,3)) ->
    os.res
   os.res
}




os.res <- datClean(NSCLC_fpkmanno,NSCLC_clin_tab,mygene)

os.res %>% 
  filter(gene%in%gene_new) %>% 
  pull(gene)->
  gene_newos

os.res %>% 
  filter(HR>1)->
  os_res_danger

os.res %>% 
  filter(HR<1)->
  os_res_protect


# DFS time ----------------------------------------------------------------
DFS_dat <- data.table::fread("dat_input/Survival_SupplementalTable_S1_20171025_xena_sp",data.table=F)


DFS_dat %>% 
  rename("patient_id"=`_PATIENT`) %>% 
  filter(`cancer type abbreviation`%in%c("LUAD","LUSC")) %>% 
  select(patient_id,DFI,DFI.time) %>% 
  drop_na() %>% 
  arrange(desc(DFI.time)) %>% 
  distinct(patient_id,.keep_all = T) %>% 
  column_to_rownames("patient_id")->
  DFS.dat



datClean_sec <- function(fpkmanno,clin_tab,mygene){
  
  mygene <- mygene[which(mygene%in%rownames(fpkmanno))] 
  
  
  fpkmanno %>% 
    t() %>% 
    as.data.frame() %>% 
    dplyr::mutate(id=substr(rownames(.),1,12)) %>% 
    dplyr::mutate(ind=rowSums(across(where(is.numeric)))) %>% 
    dplyr::arrange(desc(ind)) %>% 
    dplyr::distinct(id,.keep_all = T) %>% 
    dplyr::select(id,any_of(mygene)) %>% 
    dplyr::inner_join(clin_tab %>% 
                        rownames_to_column("id")) %>% 
    column_to_rownames("id") %>% 
    dplyr::select(any_of(sort(mygene)),everything()) %>% 
    dplyr::rename("surv_time"=DFI.time,"surv_status"=DFI) %>% 
    dplyr::select(surv_time,surv_status,everything())->
    dat_clean 
  
  
  map_df(names(dat_clean)[3:which(names(dat_clean)==sort(mygene)[length(mygene)])],
         simple_sur,
         dat_clean) %>% 
    filter(`Log-rank p`<0.05&`Cox p`<0.05) %>% 
    mutate(across(-c(gene,high_pos,low_pos),as.numeric)) %>% 
    mutate(across(where(is.numeric),round,3)) ->
    dfs.res
}

DFS.dat$DFI
DFS.dat$DFI.time
datClean_sec(NSCLC_fpkmanno,DFS.dat,motor_gene)

datClean_sec(NSCLC_fpkmanno,DFS.dat,mtor_hubgene)


dfs.res <- datClean_sec(NSCLC_fpkmanno,DFS.dat,mygene)

dfs.res %>% 
  filter(HR>1)->
  dfs_res_danger

dfs.res %>% 
  filter(HR<1)->
  dfs_res_protect




c(dfs.res %>% 
    filter(HR>1) %>% 
    pull(gene),
  os.res %>% 
    filter(HR>1) %>% 
    pull(gene)) %>% n_distinct()



setwd("dat_output")
save(dfs.res,os.res,file="dat_output/pi3k-akt-mtor pathway gene survival analysis results.Rdata")
