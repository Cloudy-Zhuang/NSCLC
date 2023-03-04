rm(list = ls())
library(tidyverse)
library(survival)
library(survminer)

# bath uni-cox ------------------------------------------------------------
## read data

load("dat_input/NSCLC-exp-clin-data.Rdata")
load("dat_output/Mtor通路研究基因.Rdata")
## Interest gene

interest_gene <- sort(unique(c(mtor_hubgene,motor_gene)))
interest_gene <- sort(replace(interest_gene,interest_gene=="GRK2","ADRBK1"))
gene_publish[!gene_publish%in%interest_gene]
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

##9-22 [ADD GENE] FRAP1 SPP1、EFNA2、ITGA6、ITGB8、ITGB4、IL2RB、TNC、KLF5、HRH3
gene_new <- c("MTOR" , "SPP1" , "EFNA2", "ITGA6" ,"ITGB8",
              "ITGB4" ,"IL2RB" ,"TNC"  , "HRH3" )


gene_publish <- c(gene_publish,gene_new)

gene_publish <- sort(gene_publish)




gene.all <- sort(unique(c(interest_gene,gene_publish)))
#pathway.gene.all <- gene.all
#save(pathway.gene.all,file = "dat_output/mtor Clean version of pathway gene.Rdata" )
# downgene <- c("GREM1","CCL17","TSLP","IL9","CXCL12","IL15")#6个
# upgene <- c("CTF1","PGF","FIGF","FGF2")

mygene <- gene.all 
# create the NSCLC log2(fpkm+1) expression matrix
NSCLC_fpkmanno <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)
NSCLC_clin_tab <- rbind(LUAD_clin_tab,LUSC_clin_tab)

## import function
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


# up gene
# datClean(NSCLC_fpkmanno,NSCLC_clin_tab,upgene)
# datClean(LUAD_fpkmanno,LUAD_clin_tab,upgene)
# datClean(LUSC_fpkmanno,LUSC_clin_tab,upgene)
# 
# datClean(NSCLC_fpkmanno,NSCLC_clin_tab,downgene)
# datClean(LUAD_fpkmanno,LUAD_clin_tab,downgene)
# datClean(LUSC_fpkmanno,LUSC_clin_tab,downgene)


# datClean(NSCLC_fpkmanno,NSCLC_clin_tab,motor_gene)
# 
# datClean(NSCLC_fpkmanno,NSCLC_clin_tab,mtor_hubgene)

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

list("OS_res"=os.res,
     "DFS_res"=dfs.res) %>% 
  openxlsx::write.xlsx("附件3修改.xlsx")


c(dfs.res %>% 
    filter(HR>1) %>% 
    pull(gene),
  os.res %>% 
    filter(HR>1) %>% 
    pull(gene)) %>% n_distinct()



setwd("D:/桌面重要文件夹/论文写作/文章写作预分析/dat_output")
save(dfs.res,os.res,file="pi3k-akt-mtor通路基因生存分析结果.Rdata")


#Save the original calculation matrix----------------------------------------------------------------
#Save (dfs. data, os. data, file="dat_output/random forest preparation data. Rdata")
#



load("dat_output/pi3k-akt-mtor通路基因生存分析结果.Rdata")
gene.important <- c("ADCY2", "ARHGAP11A", "CAMK4", "CCNA2", "CCNB2", "CDCA5", "CENPH", "CFL1", "DSCC1", "FASLG", "MAD2L1", "MAPK10", "MKNK1", "PIK3C3", "PIK3R1", "PIK3R3", "PLK4", "PRKCB", "RAC1", "RAF1", "RPS6KA3", "TPX2")

os.res %>% 
  filter(gene%in%gene.important)->ts.os

dfs.res %>% 
  filter(gene%in%gene.important)->ts.dfs


ts.os %>% 
  mutate(index="OS") %>% 
  bind_rows(ts.dfs %>% 
    mutate(index="DFS")) %>% 
  select(index,everything()) %>% 
  select(-high_pos,-low_pos) %>% 
  select(-`Cox p`) %>% 
  mutate(`95%CI`=paste0(`Lower95%`,"-",`Upper95%`)) %>% 
  select(-c(`Lower95%`,`Upper95%`,coef)) %>% 
  openxlsx::write.xlsx("table/22基因生存分析结果表格.xlsx")

if (F) {
  gene.tab1 <- openxlsx::read.xlsx("图表整理/附件1 研究的基因列表.xlsx",sheet=1)
  gene.tab2 <- openxlsx::read.xlsx("图表整理/附件1 研究的基因列表.xlsx",sheet=2)
  gene.tab3 <- openxlsx::read.xlsx("图表整理/附件1 研究的基因列表.xlsx",sheet=3)
  gene.tab1$Description <-  str_remove_all(gene.tab1$Description,"\\[.*\\]")
  gene.tab2$Description <-  str_remove_all(gene.tab2$Description,"\\[.*\\]")
  gene.tab3$Description <-  str_remove_all(gene.tab3$Description,"\\[.*\\]")
  
 gene_input <- c("ADCY2", "ARHGAP11A", "CAMK4", "CCNA2", "CCNB2", "CDCA5", "CENPH", "CFL1", "DSCC1", "FASLG", "MAD2L1", "MKNK1", "PIK3R1", "PIK3R3", "PLK4", "PRKCB", "RAC1", "RAF1", "TPX2", "ADCY2", "MAPK10", "PIK3C3", "PIK3R3", "RAC1")
  
 gene.tab1 %>% 
   select(1,8) %>% 
   rename("gene"=1) %>% 
   rbind( gene.tab2 %>% 
            select(1,8) %>% 
            rename("gene"=1)) %>% 
   rbind( gene.tab3 %>% 
            select(1,8) %>% 
            rename("gene"=1)) %>% 
   filter(gene%in%gene_input) %>% 
   openxlsx::write.xlsx("图表整理/gene_ref.xlsx")
  list(gene.tab1,gene.tab2,gene.tab3) %>% 
    openxlsx::write.xlsx("图表整理/附件1 研究的基因列表(更新).xlsx")

   
}


library(tidyverse)

dfs.res %>% 
  filter(HR>1) %>% 
  pull(gene) %>% 
  n_distinct()
#27 
dfs.res %>% 
  filter(HR<1) %>% 
  pull(gene) %>% 
  n_distinct()
#50

os.res %>% 
  filter(HR<1) %>% 
  pull(gene) %>% 
  n_distinct()
#21
os.res %>% 
  filter(HR>1) %>% 
  pull(gene) %>% 
  n_distinct()
#120

## 交集
dfs.res %>% 
  filter(HR>1) %>% 
  pull(gene) %>% 
  intersect(os.res %>% 
              filter(HR>1) %>% 
              pull(gene)) %>% 
  n_distinct()


dfs.res %>% 
  filter(HR<1) %>% 
  pull(gene) %>% 
  intersect(os.res %>% 
              filter(HR<1) %>% 
              pull(gene)) %>% 
  n_distinct()




# Single factor forest map ------------------------------------------------------------------
load("dat_output/通路基因生存分析及表达谱分析结果.Rdata")
gene.important_sel <- setdiff(gene.important,"RPS6KA3") %>% 
  sort()
# Generate data
Creat_osdat<- function(fpkmanno,clin_tab,mygene){
  
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
  dat_clean 
}

Creat_dfsdat <- function(fpkmanno,clin_tab,mygene){
  
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
  dat_clean 
}

os.dat <- Creat_osdat(NSCLC_fpkmanno,NSCLC_clin_tab,gene.important_sel) %>% 
  select(1:2,gene.important_sel)
            
dfs.dat <- Creat_dfsdat(NSCLC_fpkmanno,DFS.dat,gene.important_sel)
  

dfs.res.sel <- filter(dfs.res,gene%in%names(dfs.dat)[3:ncol(dfs.dat)])
os.res.sel <- filter(os.res,gene%in%names(os.dat)[3:ncol(os.dat)])

os.dat.sel <- select(os.dat,1:2,all_of(os.res.sel$gene))
dfs.dat.sel <- select(dfs.dat,1:2,all_of(dfs.res.sel$gene))

### Draw single-factor forest map
if(T){
dat <- dfs.dat.sel

cox_res <-dfs.res.sel
i=3
draw_forest <- function(dat, 
                        cox_res,
                        method=c("single","multi")){

  ### Build groups based on the best truncation value
  mylist <- list()

  
  for (i in 3:ncol(dat)) {
  
    # Fixed value
    cohort <- seq(1,sum(cox_res$Num_Low[1],cox_res$Num_High[1]))
    # Retrieve the grouping in the original function
    high_pos <- as.numeric(unlist(strsplit(cox_res$high_pos[i-2],"_")))
    group <- ifelse(cohort%in% high_pos,"high","low")
    
    mylist[[i-2]] <-group 
  }
  
  # Remove NULL
  ### Generate grouping data
  mylist %>% 
    purrr::reduce(cbind) %>% 
    `colnames<-`(names(dat)[3:ncol(dat)]) %>% 
    as.data.frame() %>% 
    mutate_all(factor,level=c("low","high")) %>% 
    cbind(select(dat,1:2),.)->
    dat_group
  
  
  
  require(ezcox)

  method_set <- match.arg(method)
  method_set <- switch ( method_set ,
                         "single" = "single",
                         "multi"="multi")
  
  if ( method_set=="single") {
    
    

    mds <- get_models(mymodel_sig)
    #str(mds, max.level = 1)
    show_models(mds,merge_models = TRUE, 
                drop_controls = T)
    
  
   
  }else{
    
    ### Draw multi-factor forest map
    mymodel_mul <- ezcox(dat,
                         time="surv_time",
                         status = "surv_status",
                         covariates = names(dat_group)[4:ncol(dat_group)],
                         controls = names(dat_group)[3],
                         return_models = TRUE,
                         verbose=F)
     
    mds <- get_models( mymodel_mul)
    #str(mds, max.level = 1)
    show_models(mds,
                merge_models = TRUE, 
                drop_controls =F)
    
   
   # mymodel_mul
  }
}
}

draw_forest(dfs.dat.sel,dfs.res.sel,method = "multi")

draw_forest(os.dat.sel,os.res.sel,method = "single")
x11()

draw_forest 
  
  ## Or a function
  if (F) {
    show_forest(dat_group, 
                time="surv_time",
                status = "surv_status",
                covariates = names(dat_group)[3],
                controls =names(dat_group)[4:9],
                merge_models = T,
                drop_controls = F)
  }
  
  ## another method: ggforest
  if (F) {
    
    cause_mod <- paste(names(dat_group)[3:9],collapse = "+")
    my_formula <- as.formula(paste('Surv(surv_time, surv_status)~',"+", cause_mod))
    ### Multi-factor forest map
    # ggforest(coxph(my_model,dat_group), data =dat_group, 
    #          main = "Hazard ratio", 
    #          cpositions = c(0.10, 0.22, 0.4), 
    #          fontsize = 1.2, 
    #          refLabel = "1",
    #          noDigits = 4)
    
    
    my_model <- coxph(my_formula,dat_group)
    
    
  }
  

## Ring chart
### Genes used in mapping
posgene <- c("ARHGAP11A", "CCNA2", "CCNB2", "CDCA5", "CENPH", "CFL1", "DSCC1", "FASLG", "MAD2L1", "PLK4", "RAC1", "TPX2")
neggene <- c("ADCY2", "MAPK10", "PIK3R1", "PIK3R3", "RAF1")


load("dat_output/通路基因生存分析及表达谱分析结果.Rdata")
gene.important.tab %>% 
  filter(index=="OS") %>% 
  filter(HR>1&p<0.05) %>% 
  distinct(gene,.keep_all = T) %>% 
  select(gene,HR,`Cox p`) %>% 
  rename("p"=3) %>% 
  mutate(OS_risk=1) %>% 
  mutate(OS_pro=0) %>% 
  mutate(DFS_risk=0) %>% 
  mutate(DFS_pro=0)->OS.riskgene
# RAC1 is harmful to both DFS and OS
OS.riskgene[OS.riskgene$gene=="RAC1","DFS_risk"]=1

gene.important.tab %>% 
  filter(index=="OS") %>% 
  filter(HR<1&p<0.05) %>% 
  distinct(gene,.keep_all = T) %>% 
  select(gene,HR,`Cox p`) %>% 
  rename("p"=3) %>% 
  mutate(OS_risk=0) %>% 
  mutate(OS_pro=1) %>% 
  mutate(DFS_risk=0) %>% 
  mutate(DFS_pro=0)->OS.progene

# "ADCY2" and "PIK3R3" protect both
OS.progene[OS.progene$gene== "ADCY2" ,"DFS_pro"]=1
OS.progene[OS.progene$gene== "PIK3R3" ,"DFS_pro"]=1


gene.important.tab %>% 
  filter(index=="DFS") %>% 
  filter(HR>1&p<0.05) %>% 
  distinct(gene,.keep_all = T) %>% 
  select(gene,HR,`Cox p`) %>% 
  rename("p"=3) %>% 
  mutate(OS_risk=0) %>% 
  mutate(OS_pro=0) %>% 
  mutate(DFS_risk=1) %>% 
  mutate(DFS_pro=0)->DFS.riskgene

DFS.riskgene[DFS.riskgene$gene=="RAC1","OS_risk"]=1



gene.important.tab %>% 
  filter(index=="DFS") %>% 
  filter(HR<1&p<0.05) %>% 
  distinct(gene,.keep_all = T) %>% 
  select(gene,HR,`Cox p`) %>% 
  rename("p"=3) %>% 
  mutate(OS_risk=0) %>% 
  mutate(OS_pro=0) %>% 
  mutate(DFS_risk=0) %>% 
  mutate(DFS_pro=1)->DFS.progene

DFS.progene[DFS.progene$gene== "ADCY2" ,"OS_pro"]=1
DFS.progene[DFS.progene$gene== "PIK3R3" ,"OS_pro"]=1

dat_ggcor <- rbind(OS.progene,
                   OS.riskgene,
                   DFS.progene,
                   DFS.riskgene)

dat_ggcor %>% 
  distinct(gene,.keep_all = T) %>% 
  filter(gene!="RPS6KA3") %>% 
  mutate(gene_all=1)->
  dat_ggcor_fin


data_bind <- dat_ggcor_fin



if (F) {
  
  head(data_bind)
  #The coordinates of the three rings are slightly higher than the maximum value of log2FC
  h <- max(data_bind$HR) + 2
  tile_coord <- data.frame(gene = data_bind$gene, 
                           x = 1:nrow(data_bind), # X is the middle x coordinate
                           y1 = h, y2 = h + 2, y3 = h + 4,# 
                           # y4 = h + 5,        #Y coordinates of the three outermost circles
                           stringsAsFactors = FALSE)
  
  Angle_space <- 8 # Set angle gap size
  Angle_just <- 90 # Angle Offset
  
  #Increase angle variable
  ##Unit angle
  unit_angle <- 360 / (nrow(tile_coord) + Angle_space + 0.5) 
  text_angle <- Angle_just - cumsum(rep(unit_angle, nrow(tile_coord)))
  tile_coord$hjust <- ifelse(text_angle < -90, 1, 0)
  tile_coord$text_angle <- ifelse(text_angle < -90, text_angle + 180, text_angle)
  
  data_result <- left_join(data_bind, tile_coord, by = "gene")
  
  p1 <-  ggplot(data_result) + 
    #First circle, DFS
    geom_tile(data = filter(data_result, DFS_risk == 1),
              aes(x = x, y = y1, width = 1, height = 2), 
              fill = "grey20", color = "lightgrey", size = .3) + 
    geom_tile(data = filter(data_result,DFS_risk == 0),
              aes(x = x, y = y1, width = 1, height = 2),
              fill = "grey96", color = "lightgrey", size = .3) + 
    
    #Second lap, OS
    geom_tile(data = filter(data_result, OS_risk == 1),
              aes(x = x, y = y2, width = 1, height = 2), 
              fill = "grey20", color = "lightgrey", size = .3) + 
    #If there is a value of 0, remove the following three lines#
    geom_tile(data = filter(data_result, OS_risk == 0),
              aes(x = x, y = y2, width = 1, height = 2),
              fill = "grey96", color = "lightgrey", size = .3) +
    
    #The third lap, GENEall
    geom_tile(data = filter(data_result, gene_all == 1),
              aes(x = x, y = y3, width = 1, height = 2), 
              fill = "grey20", color = "lightgrey", size = .3)+  
    #If there is a value of 0, remove the following three lines#
    #geom_tile(data = filter(data_result, GBM == 0),
    #          aes(x = x, y = y3, width = 1, height = 2),
    #          fill = "grey96", color = "lightgrey", size = .3) + 
    #If you need to draw more circles, continue to write
    
    #Write the gene name. The gene name with survival of 1 is black, otherwise it is light gray
    geom_text(data = filter(data_result, OS_risk == 1|DFS_risk==1),
              aes(x = x, y = y3 + 1.6, label = gene,
                  angle = text_angle, hjust = hjust), color = "black", size = 4) + 
    geom_text(data = filter(data_result, OS_pro == 1|DFS_pro == 1),
              aes(x = x, y = y3 + 1.6, label = gene,
                  angle = text_angle, hjust = hjust), color = "grey43", size = 4) + 
    
    coord_polar(theta = "x", start = 0, direction = 1) + 
    ylim(-5,13) + 
    xlim(0, nrow(data_bind) + Angle_space) +
    theme_void() 
  p1
  
  p2 <- p1 + geom_col(aes(x = x, y = HR*1.2), fill = "grey85", color = "grey80",width = 0.3)
  p2
  
  #Add the coordinate axis and scale text of log2FC
  p3 <- p2 + 
    geom_segment(x = 0, y = 0, xend = 0, yend = 2, color = "black", size = .5) + 
    geom_text(x = -.5, y = 0, label = "0_", vjust = -.3, size = 1, color = "black") + 
    geom_text(x = -.5, y = 0.6, label = "0.5_", vjust = -.3, size =0.9,color = "black") +
    geom_text(x = -.5, y = 1.2, label = "1.0_", vjust = -.3, size = 1,color = "black") +
    geom_text(x = -.5, y = 1.8, label = "1.5_", vjust = -.3, size = 1, color = "black")
  p3
  #Draw P_ Value heat map
  p4 <- p3 + geom_tile(aes(x = x, y = -1.2, 
                           width = 1, height = 2, fill = p), 
                       color = "white") + 
    scale_fill_gradient(low = "firebrick3", high = "yellow") 
  p4
  #insert text
  p5 <- p4 + geom_text(x = -2.1, y = h + 4.5, label = "All genes", size = 3, color = "black") + 
    geom_text(x = -1.1, y = h + 2.5, label = "OS Impact", size = 3, color = "black") + 
    geom_text(x = -1.8, y = h + .5, label = "DFS Impact", size = 2.5, color = "black") + 
    geom_text(x = -2.1, y = h - 2.4, label = "HR", 
              angle = 90, size = 2.2, color = "black") + 
    # geom_text(x = -3.3, y = h - 5, label = "p-value", 
    #           angle = 90, size = 1.5, color = "black") + 
    geom_text(x = 0, y = -5, label = "P value", size = 2.5, color = "black")
  p5
}
