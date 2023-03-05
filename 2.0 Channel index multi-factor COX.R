## Multivariate cox analysis
# Prepare data: NSCLC clinical data set index_clean data

library(tidyverse)
library(GSVA)
library(ggpubr)
library(survival)
library(survminer)
load("dat_input/NSCLC-exp-clin-data.Rdata")
load("dat_output/TCGA队列通路指数计算结果.Rdata")

# data pre ----------------------------------------------------------------

NSCLC_fpkmanno <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)
NSCLC_clin_tab <- rbind(LUAD_clin_tab,LUSC_clin_tab)
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


# FUNCTION ----------------------------------------------------------------
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




Creat_dat.os <- function(fpkmanno,
                         clin_tab,
                         mygene,
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


Creat_dat.os(fpkmanno =index_exp,
                       mygene = "index_minus",
                       NSCLC_clin_tab)->
  dat.os.multicox

table(dat.os.multicox$index_minus>0.108) 



 Creat_dat.dfs(fpkmanno =index_exp,
                         mygene = "index_minus",
                         DFS.dat) %>% 
  rownames_to_column("id") %>% 
  inner_join(dat.os.multicox %>% 
      select(-c(1:3)) %>% 
      rownames_to_column("id" )) %>% 
  column_to_rownames("id") ->
  dat.dfs.multicox
  



### Draw a single factor forest map
if(F){
  library(ezcox)
  dat_group <-select(dat_group,1,2,sort(names(dat_group)[3:ncol(dat_group)]))
  mymodel <- ezcox(dat_group,
                   time="surv_time",
                   status = "surv_status",
                   covariates = names(dat_group)[3:ncol(dat_group)],
                   controls = NULL,
                   return_models = TRUE)
  mds <- get_models(mymodel)
  str(mds, max.level = 1)
  show_models(mds,merge_models = TRUE, drop_controls = T)
  
  x11()
  
  ### Multifactor forest mapping
  mymodel <- ezcox(dat_group,
                   time="surv_time",
                   status = "surv_status",
                   covariates = names(dat_group)[3],
                   controls = names(dat_group)[4:ncol(dat_group)],
                   return_models = TRUE)
  mds <- get_models(mymodel)
  str(mds, max.level = 1)
  show_models(mds,merge_models = TRUE, drop_controls = F)
  
  x11()
  
  
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
    ### Multifactor forest map
    # ggforest(coxph(my_model,dat_group), data =dat_group, 
    #          main = "Hazard ratio", 
    #          cpositions = c(0.10, 0.22, 0.4), 
    #          fontsize = 1.2, 
    #          refLabel = "1",
    #          noDigits = 4)
    
    
    my_model <- coxph(my_formula,dat_group)
    
    
  }
  
  
  # # Calculate the riskscore for each patient -----------------------------------------------
  ### important reference :https://mp.weixin.qq.com/s/Fg-Nefyz8BIzN17oJ23wxQ
  #Replacement of gene grouping data: Complete grouping data is constructed
  
  
  dat_clean %>% 
    select(-1) %>%
    cbind(dat_group[3:ncol(dat_group)]) %>% 
    select(-neoadjuvant_treatment )->mydat
  
  cause_mod <- paste(names(dat_group)[3:ncol(dat_group)],collapse = "+")
  my_formula <- as.formula(paste('Surv(surv_time, surv_status)~',"+", cause_mod))
  my_model <- coxph(my_formula,dat_group)
  
  
  
  my_model %>% 
    broom::tidy() %>% 
    mutate(risk=exp(estimate)) %>%  #Extraction risk factor
    arrange(term) %>% 
    mutate(gene=str_sub(term,1,-5)) %>% 
    mutate(gene=paste0(gene," high*",round(risk,2))) %>% 
    pull(gene) %>% 
    paste(collapse = "+")
  
  riskscore=predict(my_model,type="risk",newdata=mydat)
  
  
  mydat %>% 
    select(-all_of(names(dat_group)[3:ncol(dat_group)])) %>%  #Delete the gene list
    mutate(riskscore=riskscore) %>% 
    select(surv_time,surv_status,riskscore,everything())->dat_multcox
  
  
}

# Multi-factor cox in action---------------------------------------------------------------
### riskscore is used as a continuity variable
##Delete columns with too many missing values

### Do not miss more than 60% of the data
multicox.tab <- function(dat_multcox,
                         riskscore,
                         percent=0.6,
                         cutoff=0.108){
  
  name_ind=which(names(dat_multcox)==riskscore)
  
  dat_multcox %>% 
    dplyr::rename("riskscore"=name_ind) %>% 
    na_if("")->dat_multcox
    
  
 
  # Find those missing values
  col_ind <- apply(dat_multcox,2,
                   function(x){sum(is.na(x))>(nrow(dat_multcox)*percent)}) 

  mycol <- which(col_ind)
  
  test <-  dat_multcox %>% 
    select(-all_of(mycol)) %>%   
    mutate(riskscore=ifelse(riskscore>cutoff,"High","Low")) %>% 
    mutate(riskscore=factor(riskscore,levels = c("Low","High"))) %>% 
    mutate(radiation_therapy=factor(radiation_therapy,levels = c("NO","YES",NA))) %>% 
    select(-pack_years_smoked.exposures)  
  
  
  
    pacman::p_load(gtsummary, gt, finalfit, survival,data.table)
    
    do_cox = function(test){
      
      
      data_for_cox=test
      
      # Define functions for data collation
      data_clean = function(data){
        data %>% 
          filter(str_detect(explanatory, "Unknown", negate = T)) %>% 
          mutate_at(c("HR", "L95", "U95"), ~round(.,2)) %>% # Keep two significant digits
          mutate(HR = ifelse(is.na(HR), "-", paste0(HR, "(", L95, ",", U95, ")")), # Integrate HR and 95% confidence intervals
                 p = ifelse(is.na(p), "-", p_tidy(p, 3, prefix = ""))) %>% # P-value formatting
          dplyr::select(explanatory, HR, p) %>% # 3列
          dplyr::rename("HR (95% CI)" = "HR")
      }
      
      #### COX univariate analysis ----
      # https://finalfit.org/reference/coxphuni.html
      explanatory_uni = setdiff(colnames(data_for_cox), c("surv_status", "surv_time"))
      dependent = "Surv(surv_time, surv_status)"
      coxphuni_res = data_for_cox %>%
        coxphuni(dependent, explanatory_uni) %>%
        fit2df(condense = F) %>% # 结果不合并
        data_clean() 
      
      # significant explanatory for coxphmulti
      explanatory_multi = coxphuni_res %>% 
        filter(p<0.05) %>% 
        pull(explanatory) %>% 
        str_extract(., paste0(explanatory_uni, collapse = "|")) %>% 
        unique()
      
      
      #### Result presentation ----
      # The table should correspond to the clinical variable name displayed (the length of labs should be consistent with the number of coxphuni_res rows!).
      # labs =c("riskscore(High vs.Low)",
      #         "Age(≥65 vs.<65 )",
      #         "Treatment Outcome(SD vs.PD)",
      #         "Treatment Outcome(PR vs.PD)",
      #         "Treatment Outcome(CR vs.PD)",
      #         "Cigeret packs somked(≥40/year vs. <40/year)",
      #         "M Stage (M1 vs. M0)",
      #         "N Stage (N2-N3 vs. N0-N1)",
      #         "T Stage (T3-T4 vs. T1-T2)",
      #         "Accepted Radiation Therapy (Yes vs.No )",
      #         "Gender (Female vs. Male)",
      #         "Clinical Stage (III-IV vs. I-II)",
      #         "Smoking exposure (≥20 years vs <20 years)"
      #        )
      labs=coxphuni_res$explanatory
      #### COX multivariate analysis ----
      # https://finalfit.org/reference/coxphmulti.html
      coxphmulti_res = data_for_cox %>%
        coxphmulti(dependent, explanatory_multi) %>%
        fit2df(condense = F) %>% # Result unmerge
        data_clean()
      
      cox_res = full_join(coxphuni_res, coxphmulti_res,
                          by = "explanatory") %>% 
        mutate(explanatory = labs) %>% 
        dplyr::rename("Variables" = "explanatory") %>% 
        expss::if_na(.,"-")
      
      
      
      
      # Make a table
      # https://gt.rstudio.com/index.html
      cox_res %>% gt() %>% 
        tab_spanner(label = md("**Univariate analysis**"), columns = ends_with(".x")) %>%
        tab_spanner(label = md("**Multivariate analysis**"), columns = ends_with(".y")) %>%
        # Change the column name displayed
        cols_label("HR (95% CI).x" = "HR (95% CI)", "p.x" = md("***p***"), 
                   "HR (95% CI).y" = "HR (95% CI)", "p.y" = md("***p***")) %>%
        tab_style(style = list(cell_fill(color = "white"), # White background
                               cell_borders(color = "white")), # The table area has no border
                  locations = cells_body()) %>%
        tab_style(style = cell_text(weight = "bold"), # Prominent P-value font bold
                  locations = list(cells_body(columns = c(p.x), rows = p.x < 0.05 & p.x != "-"), # If "Error: The object to 'data' is not a 'gt_tbl' object", prefix the vars function with the name gt::vars
                                   cells_body(columns = c(p.y), rows = p.y < 0.05 & p.y != "-"))) #%>%
      #gtsave(paste0(coxPath, "/", surv_type, "_cox.html")) # output
      
    }
    
    do_cox(test) 
    

  
}


multicox.tab(dat_multcox = dat.os.multicox,
             riskscore = "index_minus",
             percent = 0.5,
             cutoff = 0.108)


multicox.tab(dat_multcox = dat.dfs.multicox,
             riskscore = "index_minus",
             percent = 0.6,
             cutoff = 0.108 )

if (F) {
  #Survival outcome map

  dat.os.multicox$surv_status
  dat.os.multicox %>%
    ggplot(aes(factor(surv_status),
               index_minus,
               fill=factor(surv_status)))+
    geom_boxplot()+
    scale_fill_discrete(labels=c("Alive","Dead"),name="Status",
                        type =  c("#00A1D5FF","#B24745FF"))+
    scale_x_discrete(labels=c("Alive","Dead"))+
    theme_bw(base_size = 18)+
    theme(legend.position = "top")+
    ggpubr::stat_compare_means(size=5,
                               comparisons = list(c("1","0")),
                               label.x = 1.2,
                               label="p.signif",
                               show.legend = F,
                          symnum.args =list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                        symbols = c("***", "***", "**", "*", "ns")))+
    labs(x="",y="Pathway index")->a

  
  #Treatment outcome box diagram
  dat.os.multicox$Treatment_outcome %>% 
    na.omit() %>% 
    unique() %>% 
    as.character() %>% 
    combn(2,simplify = F)->com
  
  
  dat.os.multicox %>% 
    filter(!is.na(Treatment_outcome)) %>% 
    ggplot(aes(Treatment_outcome,
               index_minus,
               fill=Treatment_outcome))+
    geom_boxplot(alpha=0.2,width=0.3,
                 position=position_dodge(width=0.8),
                 size=0.75,outlier.colour = NA,
                 )+
    geom_violin(alpha=0.2,width=0.9,
                position=position_dodge(width=0.8),
                size=0.75
                )+
    ggpubr::stat_compare_means(comparisons = com,
                               label = "p.signif",
                               size=4)+
    ggpubr::stat_compare_means(size=5,
                               label.y = 1.8,
                               show.legend = F)+
    theme_bw(base_size = 18)+
    theme(legend.position = "top")+
    #ggsci::scale_fill_jco()+
    labs(x=" ",fill="Treatment Outcome")+
    labs(x="",y="Pathway index")->b
library(patchwork)
a+b+plot_layout(ncol=2,widths=c(1,1.5))+x11()

}
dat.os.multicox$surv_time %>% summary()
##Draw a multifactor table






#Line diagram model
if (F) {
  library(pacman)
  p_load(survival, rms, tidyverse)
  dat <- dat.os.multicox %>% 
    dplyr::select(1:2,index_minus,Treatment_outcome) %>% 
     mutate(index_minus=ifelse(index_minus>0.108,"high","low")) %>% 
     mutate(index_minus=factor(index_minus,levels = c("low","high")))
  
  myformula <- as.formula(paste0("Surv(surv_time, surv_status)~",
                                 paste(colnames(dat)[-c(1:2)],collapse = "+")))
  
  
  
  
  dd<-datadist(dat)
  options(datadist="dd")
  options(na.action="na.delete")
  coxpbc<-cph(formula =myformula ,
              data=dat,x=T,y=T,surv = T,
              na.action=na.delete,
              singular.ok =F)  #,time.inc =2920
  
  #print(coxpbc)
  surv<-Survival(coxpbc) 
  surv_1year<-function(x) surv(365,x)
  surv_3year<-function(x) surv(1095,x)
  surv_5year<-function(x) surv(1825,x)
  
  x<-nomogram(coxpbc,fun = list(surv_1year,surv_3year,surv_1year),lp=T,
              funlabel = c('1-year survival Probability','3-year survival Probability',
                           '5-year survival Probability'),
              maxscale = 100,
              fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
  
  #pdf("nomogram_classical.pdf",width = 12,height = 10)
  plot(x, lplabel="Linear Predictor",
       xfrac=.35,varname.label=TRUE, 
       varname.label.sep="=", ia.space=.2, 
       tck=NA, tcl=-0.20, lmgp=0.3,
       points.label='Points',
       total.points.label='Total Points',
       total.sep.page=FALSE, 
       cap.labels=FALSE,
       cex.var = 1,cex.axis = 0.8,lwd=5,
       label.every = 1,col.grid = gray(c(0.8, 0.95)))
  
  
  
  f1<-cph(formula = myformula ,data=dat,x=T,y=T,
          surv = T,
          na.action=na.delete,time.inc = 365*1) 
  
  #Parameter m=50 means that each group of 50 samples will be calculated repeatedly
  cal1<-calibrate(f1, cmethod="KM", method="boot",
                  u=365,m=20,B=1000) 
  
  f3<-cph(formula =myformula,data=dat,x=T,y=T,surv = T,
          na.action=na.delete,time.inc = 365*3) 
  cal3<-calibrate(f3, cmethod="KM", method="boot",
                  u=1095,m=20,B=1000)
  
  f5<-cph(formula =myformula,data=dat,x=T,y=T,surv = T,
          na.action=na.delete,time.inc = 365*5) 
  cal5<-calibrate(f5, cmethod="KM",
                  method="boot",
                  u=1825,m=20,B=1000)
  
  
  
  plot(cal1,
       lwd = 2,#error bar thickness
       lty = 1,#The type of error bar ranges from 0 to 6
       errbar.col = c("#2166AC"),#Color of the error bar
       xlim = c(0,1),ylim= c(0,1),
       xlab = "Nomogram-prediced OS (%)",
       ylab = "Observed OS (%)",
       cex.lab=1.2, cex.axis=1, 
       cex.main=1.2, cex.sub=0.6) #Word size
  
  lines(cal1[,c('mean.predicted',"KM")], 
        type = 'b', #The type of line, it could be "p","b","o"
        lwd = 2, #The thickness of the wire
        pch = 16, #The shape of the point, it could be 0 minus 20

        col = c("#2166AC")) #Color of line
  
  
  plot(cal3,
       lwd = 2,
       lty = 1,
       errbar.col = c("#00BFFF"),
       xlim = c(0,1),ylim= c(0,1),
       xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
       col = c("#00BFFF"),
       cex.lab=1.2,cex.axis=1, cex.main=1.2, 
       cex.sub=0.6,add = T)
  lines(cal3[,c('mean.predicted',"KM")],
        type= 'b',
        lwd = 2,
        col = c("#00BFFF"),
        pch = 16)
  
  
  #5 years
  
  plot(cal5,
       lwd = 2,
       lty = 1,
       errbar.col = c("#B2182B"),
       xlim = c(0,1),ylim= c(0,1),
       xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
       col = c("#B2182B"),
       cex.lab=1.2,cex.axis=1, cex.main=1.2, 
       cex.sub=0.6,add = T)
  lines(cal5[,c('mean.predicted',"KM")],
        type= 'b',
        lwd = 2,
        col = c("#B2182B"),
        pch = 16)
  mtext("")
  
  abline(0,1, lwd = 2, lty = 3, col = c("#224444"))
  
  legend("topleft", #The location of the legend
         legend = c("1-year","3-year","5-year"), #Legend text
         col =c("#2166AC","#00BFFF","#B2182B"), #The color of the legend line, corresponding to the text
         lwd = 2,#The thickness of the center line in the legend
         cex = 1.2,#Legend font size
         bty = "n")#Legend borders are not displayed
  
  
  
  
}




