library(tidyverse)
library(GSVA)
library(ggpubr)
library(survival)
library(survminer)
load("dat_input/NSCLC-exp-clin-data.Rdata")


# data pre ----------------------------------------------------------------

NSCLC_fpkmanno <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)
NSCLC_clin_tab <- rbind(LUAD_clin_tab,LUSC_clin_tab)

## import function
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
    data_for_surv=select(dat,all_of(gene),1,2)
   
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
    fit = surv_fit(surv.obj, data = data_for_surv) # 用survival::survfit报错: https://github.com/kassambara/survminer/issues/403
    # log-rank p http://rpkgs.datanovia.com/survminer/reference/surv_pvalue.html
    pval = surv_pvalue(fit)$pval
    # Cox proportional hazards regression model
    cox = coxph(surv.obj, data = data_for_surv)
    cox_summary = summary(cox)
    cox_res = c(cox_summary$coefficients[,c(1:2,5)],
                cox_summary$conf.int[,c(3,4)])
    
    # Extract the result and output it
    surv_res = matrix(c(gene,opticut,as.numeric(table(data_for_surv$group)), pval, cox_res), nrow = 1)
    colnames(surv_res) = c("gene","cutoff_value","Num_Low", "Num_High", "Log-rank p", "coef", "HR", "Cox p", "Lower95%", "Upper95%")
    
    
    as_tibble(surv_res)
  }
  
}



datClean.os <- function(fpkmanno,
                        clin_tab,
                        mygene,
                        method,
                        option.cutoff=NULL){
  
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
         dat_clean,
         method=method,
         option.cutoff=option.cutoff) %>% 
    filter(`Log-rank p`<0.05&`Cox p`<0.05) %>% 
    mutate(across(-gene,as.numeric)) %>% 
    mutate(across(where(is.numeric),round,3)) ->
    os.res
  os.res
}


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



datClean.dfs <- function(fpkmanno,clin_tab,
                         mygene,
                         method,
                         option.cutoff=NULL){
  
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
         dat_clean,
         method=method,
         option.cutoff=option.cutoff) %>% 
   # filter(`Log-rank p`<0.05&`Cox p`<0.05) %>% 
    mutate(across(-gene,as.numeric)) %>% 
    mutate(across(where(is.numeric),round,3)) ->
    dfs.res
}



## Defined pathway up-regulation related genes

gene.up <- sort(c("CFL1", "RAC1", "CCNA2", "PLK4", "CDCA5", "MAD2L1", "TPX2", "CCNB2", "ARHGAP11A", "DSCC1", "CENPH", "FASLG"))

## Pathway downregulates related genes
gene.down <- sort(c("PIK3R3", "ADCY2", "PIK3R1", "MAPK10", "RAF1"))
gene.up %>% 
  as.data.frame() %>% 
  set_names("gene") %>% 
  mutate(group="up") %>% 
  rbind(gene.down %>% 
          as.data.frame() %>% 
          set_names("gene") %>% 
          mutate(group="down")) %>% 
  openxlsx::write.xlsx("dat_output/通路上下调基因.xlsx")


mylist <- list(gene.up,gene.down)
names(mylist) <- c("mtor_up","mtor_down")



# gsva analysis -----------------------------------------------------------
if (F) {
  expr <- as.matrix(NSCLC_fpkmanno)
  gsva_data <- gsva(expr,
                    mylist,
                    method = "ssgsea",
                    parallel.sz=10)
  gsva_data <- as.data.frame(gsva_data )
  
  
  gsva_data %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame()->
    index_pre
  
  index_pre %>% 
    mutate(index_plus=mtor_up+mtor_down,
           index_minus=mtor_up-mtor_down,
           index_multipy=mtor_up*mtor_down,
           index_divide=mtor_up/mtor_down) %>% 
    select(-mtor_up,-mtor_down)->
    index_clean
}


#
##save(index_clean,file = "dat_output/TCGACalculation result of queue path index.Rdata")
load("dat_output/TCGA队列通路指数计算结果.Rdata")
if (F) {
  index_clean %>% 
    mutate(group=ifelse(substr(rownames(.),14,15)=="11","normal","cancer")) ->
    index_group
  
  index_group %>% 
    ggplot(aes(group,index_plus))+
    geom_boxplot()+
    stat_compare_means()
  
  index_group %>% 
    ggplot(aes(group,index_minus,fill=group))+
    geom_boxplot()+
    stat_compare_means(size=5,label.x.npc = 0.2)+
    ggsci::scale_fill_npg()+
    theme_bw(base_size = 18)+
    labs(y="PI3K/AKT/mTOR Pathway Index")
  x11()
  
  index_group %>% 
    ggplot(aes(group,index_multipy,fill=group))+
    geom_boxplot()+
    stat_compare_means()
  index_group %>% 
    ggplot(aes(group,index_divide))+
    geom_boxplot()+
    stat_compare_means()
}


# The effects of different indices on survival --------------------------------------------------------------

index_clean %>% 
  t() %>% 
  as.data.frame()->
  index_sur


index.os <- datClean.os(fpkmanno=index_sur,
                        clin_tab = NSCLC_clin_tab,
                        mygene =rownames(index_sur),
                        method = "best_sep")
index.os <- datClean.os(fpkmanno=index_sur,
                        clin_tab = NSCLC_clin_tab,
                        mygene =rownames(index_sur),
                        method = "median")


index_sur


index.dfs.option <- datClean.dfs(index_sur,
                                DFS.dat,
                                mygene =rownames(index_sur),
                                method = "option",
                                option.cutoff = 0.108)


index.dfs <- datClean_dfs(index_sur,
                          DFS.dat,
                          mygene =rownames(index_sur),
                          method = "best_sep")


#save(index.os,index.dfs,file = "dat_output/cox analysis results of TCGA-NSCLC index.Rdata")


# Survival analysis drawing ------------------------------------------------------------------

Creat_dat.dfs <- function(fpkmanno,clin_tab,mygene){
  
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




Creat_dat.os <- function(fpkmanno,clin_tab,mygene,
                         survival){
  
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

}





index_exp <- as.data.frame(t(index_clean))


dat.os <- Creat_dat.os(fpkmanno =index_exp,
                       mygene = c("index_minus","index_divide"),
                       NSCLC_clin_tab)






dat.dfs <- Creat_dat.dfs(fpkmanno =index_exp,
                        mygene = "index_minus",
                       DFS.dat)



# save data ---------------------------------------------------------------
save(index_exp,file="dat_output/通路指数数据.Rdata")



draw_kmplot <- function(gene_input,
                        cox_res,
                        survival=NULL,
                        dat){
  
  maxX = ceiling(max(dat$surv_time))*1.05
  
  filter(cox_res,gene%in%gene_input) %>% 
    pull(cutoff_value) %>% 
    as.numeric()->cutoff 
  
  ##cox regression and HR value extraction
  my.surv <- Surv(dat$surv_time, dat$surv_status)##K-M method survival analysis fitting time and outcome
 
  group <-ifelse(dat[[gene_input]]>cutoff,"high","low") #Set grouping
  
  group <- factor(group, levels = c("low", "high"))
  survival_dat <- data.frame(group = group)
  dat$group <-  group
  fit <- survfit(my.surv ~ group)##The fitting is done. This is fitting by group using logrank's method, and this is ready to graph
。

  
  ##Calculate HR and 95%CI
  
  ##Extraction p-value
  
  pval = filter(cox_res,gene%in%gene_input) %>% 
    pull(5) %>% 
    as.numeric() %>% 
    round(3)
  
  
  
  x = summary(coxph(Surv(surv_time, surv_status)~group, data = dat))##What cox regression fitted was the expression quantity.
  HR =filter(cox_res,gene%in%gene_input) %>% 
    pull(HR) %>% 
    as.numeric() %>% 
    round(2)
  
  
  up95 = filter(cox_res,gene%in%gene_input) %>% 
    pull(`Upper95%`) %>% 
    as.numeric() %>% 
    round(3)
  
  low95 = filter(cox_res,gene%in%gene_input) %>% 
    pull(`Lower95%`) %>% 
    as.numeric() %>% 
    round(3)
  
  HR <- paste("Hazard Ratio = ", HR, sep = "")
  CI <- paste("95% CI: ", paste(round(low95,3), round(up95,3), sep = " - "), sep = "")
  pText = ifelse(pval < 0.01, formatC(pval, digits = 2, format = "E"),
                 round(pval, 3))
  #Sequence the gene expression from high to low in order to extract the boundary expression

  
  p<- ggsurvplot(fit, 
                 data = dat,
                 xlab = "Time(Days)", 
                 ylab = "Survival Probability (%)",  # ylab = "Survival probability (%)"
                # palette = c("navyblue","firebrick1"), # Curve setting
                 palette="lancet",
                 legend.title = 'Category', 
                 legend.labs = c("Low","High"),
                 font.legend = 14, legend = c(0.8,0.8), # Legend setting
                 #pval = paste0("HR = ", round(surv_res[5], 3), "\nlog-rank p = ", pText), pval.size = 5, pval.coord = c(0.6, 18), # P-value setting
                 pval = F, # 用annotate设置文本
                 #parse = T,
                 #break.time.by = 5, # x axis scale interval
                 risk.table = F, risk.table.height = 0.2, risk.table.fontsize = 5, 
                 risk.table.col = "strata", tables.y.text = F, # Risk table setting
                 font.x = 16, 
                 font.y = 16, 
                 #font.xtickslab = 13, 
                 #font.ytickslab = 15,
                 fun = function(y) y*100
  )  # The ordinate scale is set to percent
  # Set "p" italics! 
  if (!is.null(survival)) {
    p$plot = p$plot + 
      annotate("text", x = 0.5, 
               y = c(50,40,30,20,10), 
               label = c(survival,
                         gene_input,
                         HR,
                         CI,
                         "log-rank"), 
               size = 5, 
               hjust = 0, vjust = 1) + # hjust sets left and right: 0-left align, 0.5-center, 1-right align; vjust sets up - down for 0-1
      annotate("text", 
               x = 0.5, y = 10, 
               label = as.character(as.expression(
                 substitute(italic(p)~":"~pText, 
                            list(pText = formatC(pText, 
                                                 format = "e", 
                                                 digits = 2))))), 
               size = 5, 
               hjust = -0.7, 
               vjust = 1, parse = T)
  }else{
    p$plot = p$plot + 
      annotate("text", x = 0.5, 
               y = c(40,30,20,10), 
               label = c(
                         gene_input,
                         HR,
                         CI,
                         "log-rank"), 
               size = 5, 
               hjust = 0, vjust = 1) + # hjust sets left and right: 0-left align, 0.5-center, 1-right align; vjust sets up - down for 0-1

      annotate("text", 
               x = 0.5, y = 10, 
               label = as.character(as.expression(
                 substitute(italic(p)~":"~pText, 
                            list(pText = formatC(pText, 
                                                 format = "e", 
                                                 digits = 2))))), 
               size = 5, 
               hjust = -0.7, 
               vjust = 1, parse = T)
  }
  
  p
}
library(survminer)
library(survival)
draw_kmplot(gene_input = "index_minus",
            dat=dat.os,
            cox_res = index.os,
            survival = "OS")
x11()
draw_kmplot(gene_input = "index_minus",
            dat=dat.dfs,
            cox_res = index.dfs.option,
            survival = "DFS")
x11()

