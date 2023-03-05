library(tidyverse)
library(GEOquery)
library(limma)
library(data.table)

normal_gse <- function(dat){
  exprSet <- exprs(dat)
  ex <- exprs(dat)
  
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
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
  exprSet
  
}
##GSE17710 cannot be commented and can only be commented manually
load("dat_input/肺腺癌或鳞癌数据集/GSE17710.Rdata")
gse17710_gene <- fread("dat_input/肺腺癌或鳞癌数据集/GSE17710_family.soft/GSE17710_family.soft",
                       check.names=T,skip=125)
GSE17710.dat_normal <- normal_gse(GSE17710.dat)
GSE17710.dat_normal %>% 
  rownames_to_column("ID") %>% 
  inner_join(gse17710_gene %>% 
             select("ID","ORF") %>% 
               mutate(ID=as.character(ID))) %>% 
  mutate(index=rowMeans(across(where(is.numeric)))) %>% 
  arrange(desc(index)) %>% 
  distinct(ORF,.keep_all = T) %>% 
  select(-ID) %>% 
  column_to_rownames("ORF") ->
  GSE17710.dat_clean
  
GSE17710.pheno %>% 
  select(2,58:73)->
  GSE17710.pheno_clean


save( GSE17710.dat_clean,GSE17710.pheno_clean,
      file = "GSE17710_clean.Rdata")



load("GSE17710_clean.Rdata")


dat_list <- list.files("dat_input/肺腺癌或鳞癌数据集/")
dat_list[2]


# GSE26939 ----------------------------------------------------------------
load("dat_input/肺腺癌或鳞癌数据集/GSE26939.Rdata")



GSE26939.dat_normal <-  normal_gse(GSE26939.dat)
  
  

gse26939_gene <- fread("dat_input/肺腺癌或鳞癌数据集/GSE26939_family.soft.gz",
                       check.names=T,skip=186)



GSE26939.dat_normal %>% 
  rownames_to_column("ID") %>% 
  inner_join(gse26939_gene  %>% 
               select("ID","ORF") %>% 
               mutate(ID=as.character(ID))) %>% 
  mutate(index=rowMeans(across(where(is.numeric)))) %>% 
  arrange(desc(index)) %>% 
  distinct(ORF,.keep_all = T) %>% 
  select(-ID) %>% 
  column_to_rownames("ORF") ->
  GSE26939.dat_clean

GSE26939.pheno %>% 
  select(2,65:85)->
  GSE26939.pheno_clean

save(GSE26939.dat_clean, GSE26939.pheno_clean,
     file =  "GSE26939_clean.Rdata")



# GSE GSE63459 ---------------------------------------------------------------------
load("dat_input/肺腺癌或鳞癌数据集/GSE63459.Rdata")

GSE63459.dat.normal <- normal_gse(GSE63459.dat)

gse63459_gene <- fread("dat_input/肺腺癌或鳞癌数据集/GSE63459_family.soft.gz",
                       check.names=T,skip=157)


GSE63459.dat.normal  %>% 
  rownames_to_column("ID") %>% 
  inner_join(gse63459_gene  %>% 
               select(1,6)) %>% 
  mutate(index=rowMeans(across(where(is.numeric)))) %>% 
  arrange(desc(index)) %>% 
  distinct(ILMN_Gene,.keep_all = T) %>% 
  select(-ID) %>% 
  column_to_rownames("ILMN_Gene") ->
  GSE63459.dat_clean

GSE63459.pheno %>% 
  select(2,47:61)->
  GSE63459.pheno_clean

save( GSE63459.dat_clean,   GSE63459.pheno_clean,
     file =  "GSE63459_clean.Rdata")


# GSE GSE63805 ---------------------------------------------------------------------
load("dat_input/肺腺癌或鳞癌数据集/GSE63805.Rdata")

GSE63805.dat_normal <- normal_gse(GSE63805.dat)
# This data is miRNA, unqualified, skip


# GSEGSE72094 ---------------------------------------------------------------------
load("dat_input/肺腺癌或鳞癌数据集/GSE72094.Rdata")

GSE72094.dat_normal <- normal_gse(GSE72094.dat)

gse72094_gene <-fread("dat_input/肺腺癌或鳞癌数据集/GSE72094_family.soft.gz",
                       check.names=T,skip=508)


gse72094_gene$GeneSymbol

GSE72094.dat_normal   %>% 
  rownames_to_column("ID") %>% 
  inner_join(gse72094_gene  %>% 
               select(1,4) %>% 
               separate_rows(GeneSymbol,sep=" ")) %>% 
  mutate(index=rowMeans(across(where(is.numeric)))) %>% 
  arrange(desc(index)) %>% 
  distinct(GeneSymbol,.keep_all = T) %>% 
  select(-ID) %>% 
  column_to_rownames("GeneSymbol") ->
  GSE72094.dat_clean

GSE72094.pheno %>% 
  select(2,49:70)->
  GSE72094.pheno_clean

save(GSE72094.dat_clean,  
      GSE72094.pheno,
      file ="GSE72094_clean.Rdata")



# GSE47115 ----------------------------------------------------------------

load("dat_input/肺腺癌或鳞癌数据集/GSE47115.Rdata")

GSE47115.dat_normal <- normal_gse(GSE47115.dat)


gse47115_gene <-fread("dat_input/肺腺癌或鳞癌数据集/GSE47115_family.soft.gz",
                      check.names=T,skip=130)



GSE47115.dat_normal   %>% 
  rownames_to_column("ID") %>% 
  inner_join(gse47115_gene   %>% 
               select(1,12)) %>% 
  mutate(index=rowMeans(across(where(is.numeric)))) %>% 
  arrange(desc(index)) %>% 
  distinct(Symbol,.keep_all = T) %>% 
  select(-ID) %>% 
  column_to_rownames("Symbol") ->
  GSE47115.dat_clean



GSE47115.pheno %>% 
  select(2,45:55)->
  GSE47115.pheno_clean


save(GSE47115.dat_clean,  
     GSE47115.pheno,
     file ="GSE47115_clean.Rdata")


# NSE quantization -------------------------------------------------------------------

## Data set quantization and survival analysis
library(tidyverse)
library(GSVA)


dat_list <- dir(path = "dat_output/肺鳞癌-腺癌手工整理数据集/",
                full.names = T,recursive=T )


for (i in seq_along(dat_list)) {
  load(dat_list[i])
}


## Defined pathway up-regulation related genes
gene.up <- sort(c("CFL1", "RAC1", "CCNA2", "PLK4", "CDCA5", "MAD2L1", "TPX2", "CCNB2", "ARHGAP11A", "DSCC1", "CENPH", "FASLG"))

## Pathway downregulates related genes
gene.down <- sort(c("PIK3R3", "ADCY2", "PIK3R1", "MAPK10", "RAF1"))

mylist <- list(gene.up,gene.down)
names(mylist) <- c("mtor_up","mtor_down")


gse.list <- ls(pattern = "dat_clean$")



for (i in seq_along(gse.list)) {
  
  
  expr <- as.matrix(get(gse.list[i]))
  
  gsva_data <- gsva(expr,
                    mylist,
                    method = "ssgsea",
                    parallel.sz=12)
  
  file.name <- str_remove(gse.list[i],".dat_clean")
  assign(paste0(file.name,".NES"),gsva_data)
  
}
save(list = ls(pattern = "NES"),
     file = "geo_手工整理单独量化结果.Rdata")









dat.list <- ls(pattern = "NES$")


for (i in seq_along(dat.list)) {
  
  
  get(dat.list[i]) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(index_minus=mtor_up-mtor_down) %>% 
    select(-mtor_up,-mtor_down)->
    index_clean
  
  assign(paste0(str_remove(dat.list[i],".NES"),
                ".gsva"),index_clean)
  
}


# Collate clinical information from GEO data ------------------------------------------------------------
ls(pattern = "pheno") %>% sort()

# GSE17710 OS

GSE17710.pheno_clean %>% 
 # select(starts_with("survival"),geo_accession,
 #        everything()) %>% 
 select(9:10,geo_accession,
         everything()) %>% 
 rename("surv_time"=1,
        "surv_status"=2) %>% 
 mutate(surv_status=as.numeric(surv_status)) %>% 
 mutate(surv_time=as.numeric(surv_time)*30) %>%     # month*30
 inner_join(GSE17710.gsva %>% 
            rownames_to_column("geo_accession")) %>% 
  column_to_rownames("geo_accession")->GSE17710_OS
  


# GSE26939 ----------------------------------------------------------------


GSE26939.pheno_clean %>% 
  select(starts_with("survival"),geo_accession,
         everything()) %>% 
  rename("surv_time"=1,
         "surv_status"=2) %>% 
  mutate(surv_time=as.numeric(surv_time)*30) %>%     # month*30
  inner_join(GSE26939.gsva %>% 
               rownames_to_column("geo_accession")) %>% 
  column_to_rownames("geo_accession")->GSE26939_OS



GSE47115.pheno %>% 
  select(2,45:55) %>% 
  select(9,5,everything()) %>% 
  rename("surv_time"=1, "surv_status"=2) %>% 
  mutate(surv_time=as.numeric(surv_time)*30) %>%# month*30
  inner_join(GSE47115.gsva %>% 
            rownames_to_column("geo_accession")) %>% 
  mutate(surv_status=ifelse(surv_status=="yes",1,0))->
  GSE26939_OS_type1


GSE47115.pheno %>% 
  select(2,45:55) %>% 
  select(9,5,everything()) %>% 
  rename("surv_time"=1, "surv_status"=2) %>% 
  mutate(surv_time=as.numeric(surv_time)*30) %>%# month*30
  inner_join(GSE47115.gsva %>% 
               rownames_to_column("geo_accession")) %>% 
  mutate(surv_status=ifelse(surv_status=="yes",0,1))->
  GSE26939_OS_type2

# GSE63459.pheno_clean ----------------------------------------------------------


GSE63459.pheno_clean %>% 
  select(14,4,everything()) %>% 
  rename("surv_time"=1, "surv_status"=2) %>% 
  mutate(surv_status=ifelse(surv_status=="yes",1,0)) %>% 
  inner_join(GSE63459.gsva %>% 
               rownames_to_column("geo_accession"))->
  GSE63459_OS



# GSE72094.pheno ----------------------------------------------------------

GSE72094.pheno %>% 
  select(2,49:70) %>% 
  select(19,23,everything()) %>% 
  rename("surv_time"=1, "surv_status"=2) %>% 
  mutate(surv_status=ifelse(surv_status=="alive",0,1)) %>% 
  inner_join(GSE72094.gsva %>% 
               rownames_to_column("geo_accession"))->
  GSE72094_OS



# Survival analysis --------------------------------------------------------------------
if (T) {
  #simple_sur
  
  simple_sur <- function(gene,
                         dat,
                         method="best_sep",
                         option.cutoff=NULL){
    
    if (!method%in%c("best_sep","median","option")) {
      stop("method should be best_sep,median or option!")
    }
    
    print(gene)
    #data=sur_dat_fin
    data_for_surv=dplyr::select(dat,all_of(gene),1,2)
    
    if (method=="best_sep") {
      # Calculate the optimal classification point
      x <-try(surv_cutpoint(data_for_surv,
                            time = "surv_time",
                            event = "surv_status",
                            variables = gene,
                            minprop = 0.1), silent=TRUE)
      if ('try-error' %in% class(x)){
        
        surv_res = matrix(c(gene,"无法分组",rep(NA,8)), nrow = 1)
        colnames(surv_res) = c("gene","cutoff_value","Num_Low", "Num_High", "Log-rank p", "coef", "HR", "Cox p", "Lower95%", "Upper95%")
        return(as_tibble(surv_res))
        
      }else{
        
        surv.cut = surv_cutpoint(data_for_surv,
                                 time = "surv_time",
                                 event = "surv_status",
                                 variables = gene,
                                 minprop = 0.1)}
      
      
      # # Optimal classification value
      opticut = surv.cut$cutpoint$cutpoint
      
      
      
      data_for_surv$group = surv_categorize(surv.cut)[,gene]
      
      data_for_surv$group=factor( data_for_surv$group,levels = c("low","high"))
   
       }else if (method=="median"){
      opticut=median(data_for_surv[[gene]])
      data_for_surv$group <- ifelse(data_for_surv[[gene]]>median(data_for_surv[[gene]]),"high","low")
      data_for_surv$group=factor( data_for_surv$group,levels = c("low","high"))
    }else if (method=="option"&!is.null(option.cutoff)){
      opticut=option.cutoff
      data_for_surv$group <- ifelse(data_for_surv[[gene]]>option.cutoff,"high","low")
      data_for_surv$group=factor( data_for_surv$group,levels = c("low","high"))
    }else stop("choose the proper method!")
    
    
    
    #Another case of not being able to group
    if (any(is.na(data_for_surv$group))){
      
      surv_res = matrix(c(gene,"无法分组",rep(NA,8)), nrow = 1)
      colnames(surv_res) = c("gene","cutoff_value","Num_Low", "Num_High", "Log-rank p", "coef", "HR", "Cox p", "Lower95%", "Upper95%")
      return(as_tibble(surv_res))}
    
    
    # Survival analysis 
    # survfit object
    surv.obj = as.formula(paste0("Surv(surv_time, surv_status)~","group"))
    # K-M survival curves
    fit = surv_fit(surv.obj, data = data_for_surv) # An error is reported with survival::survfit: https://github.com/kassambara/survminer/issues/403
    # log-rank p http://rpkgs.datanovia.com/survminer/reference/surv_pvalue.html
    pval = surv_pvalue(fit)$pval
    # Cox proportional hazards regression model
    cox = coxph(surv.obj, data = data_for_surv)
    cox_summary = summary(cox)
    cox_res = c(cox_summary$coefficients[,c(1:2,5)],
                cox_summary$conf.int[,c(3,4)])
    
    # Output the corresponding group position as well, which can be used to draw the forest map

    
    high_pos <- paste(which(data_for_surv$group%in%"high"),collapse = "_")
    low_pos <-  paste(which(data_for_surv$group%in%"low"),collapse = "_")
    # Extraction result output
    surv_res = matrix(c(gene,opticut,as.numeric(table(data_for_surv$group)), pval, cox_res,high_pos,low_pos), nrow = 1)
    colnames(surv_res) = c("gene","cutoff_value","Num_Low", "Num_High", "Log-rank p", "coef", "HR", "Cox p", "Lower95%", "Upper95%","high_pos",
                           "low_pos")
    
    
    
    as_tibble(surv_res)
  }
  
}


gselist <- ls(pattern = "_OS")

for (i in seq_along(gselist)) {
  print(i)
  mydat <- get(gselist[i])
  mydat <- filter(mydat,surv_time!="NA"|surv_time!="NA")
  mydat$surv_status <- as.numeric( mydat$surv_status)
  mydat$surv_time <- as.numeric( mydat$surv_time)
 
  cox_res <- simple_sur(gene="index_minus",
                        dat=mydat,
                        option.cutoff = 0.108,
                        method = "option")
  
  assign(paste0(str_remove(gselist[i],"_OS"),
                "_coxres"),
         cox_res )
  
}

coxlist <- ls(pattern = "_coxres")


map(coxlist,get) %>% 
  map_df(rbind)->
  cox_resall


## failed。
