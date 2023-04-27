library(tidyverse)
library(survival)
library(survminer)
# pathway index cacluated by ssGSEA and the key pathway survival-related genes
load("dat_input/GSE30219_pathwayindex.Rdata")
load("dat_output/Pathway_index_of_TCGA_NSCLC_cohort.Rdata")
load("dat_output/GSE30219.clean.Rdata")


# Collate the analysis results of cloud server--------------------------------------------------------------
dat.list <- ls(pattern = "NES$")
GSE30219.pheno$`gender:ch1` %>% table()

for (i in seq_along(dat.list)) {
  
  
  get(dat.list[i]) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(index_minus=mtor_up-mtor_down) %>% 
    dplyr::select(-mtor_up,-mtor_down)->
    index_clean
  
  assign(paste0(str_remove(dat.list[i],".NES"),
                ".gsva"),index_clean)
  
}



# Collate clinical information from GEO data ------------------------------------------------------------
ls(pattern = "pheno$") %>% sort()



##GSE30219 
GSE30219.pheno %>% 
  dplyr::select(str_which(names(.),"disease free survival|relapse|tissue:ch1")) %>% 
  dplyr::rename("surv_time"=1,
         "surv_status"=2,
         "group"=3) %>% 
  filter(!is.na(surv_time)) %>% 
  filter(surv_status!="na") %>% 
  mutate(surv_time=as.numeric(surv_time)*30)->
  GSE30219.pheno.clean.dfs


GSE30219.pheno %>% 
  dplyr::select(str_which(names(.),"follow-up|status:ch1|tissue:ch1")) %>% 
  dplyr::rename("surv_time"=1,
         "surv_status"=2,
         "group"=3) %>% 
  filter(surv_time!="NTL") %>% 
  filter(surv_time!="NTL") %>% 
  mutate(surv_status=ifelse(surv_status=="ALIVE",0,1)) %>%   
  mutate(surv_time=as.numeric(surv_time)*30)->
  GSE30219.pheno.clean.os


# Survival analysis --------------------------------------------------------------------
load("dat_output/cox analysis results of TCGA-NSCLC index.Rdata")

## Construct extraction data function and drawing function
Creat_dat <- function(fpkmanno,clin_tab,mygene){
  
  mygene <- mygene[which(mygene%in%names(fpkmanno))] 
  
  fpkmanno %>% 
    rownames_to_column("id") %>% 
    dplyr::select(id,any_of(mygene)) %>% 
    dplyr::inner_join(clin_tab %>% 
                        rownames_to_column("id")) %>% 
    column_to_rownames("id") %>% 
    dplyr::select(any_of(sort(mygene)),everything()) %>% 
    dplyr::select(surv_time,surv_status,everything())->
    dat_clean 
  
  dat_clean 
  
}




draw_kmplot.GEO<- function(gene_input,
                        cutoff,
                        survival=NULL,
                        dat,
                        palette=NULL){
  
  maxX = ceiling(max(dat$surv_time))*1.05
  ##cox regression and HR value extraction
  surv.obj <- as.formula(paste0("Surv(surv_time, surv_status)~","group"))##Fitting time and outcome of K-M survival analysis
  
  group <-ifelse(dat[[gene_input]]>cutoff,"high","low") #Set grouping
  
  group <- factor(group, levels = c("low", "high"))
  survival_dat <- data.frame(group = group)
  dat$group <-  group
  
  fit <- surv_fit(surv.obj,data = dat)
  
  # log-rank p http://rpkgs.datanovia.com/survminer/reference/surv_pvalue.html
  pval = surv_pvalue(fit)$pval
  # Cox proportional hazards regression model
  cox = coxph( surv.obj, data = dat)
  cox_summary = summary(cox)
  cox_res = c(cox_summary$coefficients[,c(1:2,5)],
              cox_summary$conf.int[,c(3,4)])
 
  # Extract the result and output it

  HR <- paste("Hazard Ratio = ",round(cox_res[2],3), sep = "")
  CI <- paste("95% CI: ", paste(round( cox_res[4],3), round(cox_res[5],3), sep = " - "), sep = "")
  pText = ifelse(pval < 0.01, formatC(pval, digits = 2, format = "E"),
                 round(pval, 3))
  #Sequence the gene expression from high to low in order to extract the boundary expression
  
  p<- ggsurvplot(fit,
                 
                 xlab = "Time(Days)", 
                 ylab = "Survival Probability (%)",  # ylab = "Survival probability (%)"
                 palette =palette, # Curve setting
                 #palette="nejm",
                 #palette = "Dark2",
                 legend.title = 'Category', 
                 legend.labs = c("Low","High"),
                 font.legend = 14, legend = c(0.8,0.8), # Legend setting
                 #pval = paste0("HR = ", round(surv_res[5], 3), "\nlog-rank p = ", pText), pval.size = 5, pval.coord = c(0.6, 18), # P-value setting
                 pval = F, # Set the text in annotate
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




#One-by-one verification -----------------------------------------------------------------
mytry <-function(dat.input,cutoff,survival){
  
  ind <- str_extract(dat.input,"(^GSE[0-9]+)")
  
  dat <- Creat_dat(get(paste0(ind,".gsva")),
                   get(paste0(ind,".pheno.clean.",survival)), 
                   "index_minus") %>% 
    select(-group) %>% 
    mutate(surv_status=as.numeric(surv_status))
  
  ##cox regression and HR value extraction
  my.surv <- Surv(dat$surv_time, dat$surv_status)##Fitting time and outcome of K-M survival analysis
  
  
  group <-ifelse(dat[["index_minus"]]>cutoff,"high","low") #Set grouping
  
  group <- factor(group, levels = c("low", "high"))
  
  
  # Survival analysis 
  # survfit object
  surv.obj = as.formula(paste0("Surv(surv_time, surv_status)~","group"))
  # K-M survival curves
  fit = surv_fit(surv.obj, data = dat) 
  # log-rank p http://rpkgs.datanovia.com/survminer/reference/surv_pvalue.html
  pval = surv_pvalue(fit)$pval
  # Cox proportional hazards regression model
  cox = coxph(surv.obj, data = dat)
  cox_summary = summary(cox)
  cox_res = c(cox_summary$coefficients[,c(1:2,5)],
              cox_summary$conf.int[,c(3,4)])
  
  # Extract the result and output it
  surv_res = matrix(c(ind ,cutoff,as.numeric(table(group)), pval, cox_res), nrow = 1)
  colnames(surv_res) = c("GSE","cutoff_value","Num_Low", "Num_High", "Log-rank p", "coef", "HR", "Cox p", "Lower95%", "Upper95%")
  
  
  as_tibble(surv_res)
  
}



# Optimum cut-off value ------------------------------------------------------------------


datlist <- ls(pattern = "clean\\.os$")
cutoff <- 0.108





map_df(datlist,mytry,cutoff,"os") %>% 
  mutate(across(-1,as.numeric)) %>%  
  mutate(across(where(is.numeric),round,3)) %>% 
  DT::datatable()




# drawing ----------------------------------------------------------------------
## When the OS median value is 0.377

dat1 <- Creat_dat(GSE30219.gsva,
                  GSE30219.pheno.clean.os,
                  "index_minus") %>% 
  select(-group)


draw_kmplot.GEO(gene_input = "index_minus",
                dat = dat1,
                cutoff = 0.108,
                survival = "OS",
                palette = c("#0a719b","#9d0f26"))->p


x11()
p


## Multifactor table supplement
# Write out clinical data for manual screening


openxlsx::write.xlsx(GSE30219.pheno,
                     "dat_output/gse30219 Clinical information.xlsx")
#Manually delete irrelevant information
gse30219.clin <- openxlsx::read.xlsx("dat_output/gse30219 Clinical information.xlsx")



which(gse30219.clin=="NTL",arr.ind = T)
gse30219.clin $`age.at.surgery:ch1` %>% table()
gse30219.clin %>% 
  mutate_all(trimws) %>% 
  mutate(across(everything(),~na_if(.,"NTL"))) %>% 
  mutate(across(everything(),~na_if(.,"ND"))) %>% 
  mutate(age=as.numeric(`age.at.surgery:ch1`)) %>% 
  mutate(age=ifelse(age>=65,"≥65","<65")) %>% 
  mutate(age=factor(age,levels =c("<65","≥65"))) %>% 
  mutate(gender=factor(`gender:ch1`)) %>% 
  mutate(histo=fct_collapse(`histology:ch1`,
                            LUAD="ADC",
                            LUSC="SQC",
                            Others=setdiff(gse30219.clin$`histology:ch1`,c("ADC","SQC"))[-1])) %>% 
  mutate(stage=case_when(`pm.stage:ch1`=="M1"~"IV",
                         `pn.stage:ch1`=="N3"&`pm.stage:ch1`!="M1"~"III",
                         `pn.stage:ch1`=="N2"&`pm.stage:ch1`!="M1"~"III",
                         `pt.stage:ch1`%in%c("T1","T2")&`pn.stage:ch1`=="N0"&`pm.stage:ch1`!="M1"~"I",
                         `pt.stage:ch1`%in%c("T1","T2")&`pn.stage:ch1`=="N1"&`pm.stage:ch1`!="M1"~"II",
                         `pt.stage:ch1`%in%c("T3","T4")&`pn.stage:ch1`%in%c("N1","N2")&`pm.stage:ch1`!="M1"~"III",
                         `pt.stage:ch1`=="T3"&`pn.stage:ch1`=="N0"&`pm.stage:ch1`!="M1"~"II")) %>% 
  mutate(stage=case_when(stage%in%c("I","II")~"I-II",
                         stage%in%c("III","IV")~"III-IV")) %>% 
  mutate(stage=factor(stage,levels = c("I-II","III-IV"))) %>% 
  mutate(status=factor(`status:ch1`)) %>% 
  dplyr::select(geo_accession,age,stage,
         gender,histo) %>% 
  mutate(histo=factor(histo,levels = c("Others","LUAD","LUSC"))) %>% 
  inner_join(GSE30219.gsva %>% 
              rownames_to_column("geo_accession")) %>% 
  inner_join(GSE30219.pheno.clean.os %>% 
               rownames_to_column("geo_accession")) %>% 
  mutate(pathway_index=ifelse(index_minus>=0.108,"High","Low")) %>%
  mutate(pathway_index=factor(pathway_index,levels = c("Low","High"))) %>% 
  dplyr::select(-group,index_minus,-geo_accession) %>% 
  dplyr::select(surv_time,surv_status,pathway_index,everything())->
  dat_group

dat_group$stage %>% 
  as.character() ->dat_group$stage
dat_group$stage[153] <- "I-II"  #T2N1M0 referring to gse30219.clin
dat_group$stage <- factor(dat_group$stage,
                          levels =  c("I-II","III-IV"))
  
which(is.na(dat_group$age))#68   referring to gse30219.clin
dat_group$age <- as.character(dat_group$age)
dat_group$age[153] <- "≥65"
dat_group$age <- factor(dat_group$age,
                    levels =c("<65","≥65"))

which(is.na(dat_group$histo))
dat_group$histo %>% table()
gse30219.clin$`histology:ch1` %>% table()
table(dat_group$pathway_index)
table(dat_group$pathway_index,dat_group$stage)

#ploting
if (F) {
  dat_group %>% 
    drop_na() %>% 
    ggplot(aes(stage,index_minus,fill=stage))+
    geom_boxplot(alpha=0.8)+
    ggpubr::stat_compare_means(comparisons = list(c("I-II","III-IV")),
                               label="p.signif")+
    theme_bw(base_size = 18)+
    theme(legend.position = "top")+
    scale_fill_manual(values = c("#ffffb3","#8dd3c7"))+
    labs(y="Pathway index")+
    x11()
  
  library(RColorBrewer)
  display.brewer.all()
  dat_group %>% 
    drop_na() %>% 
    mutate(status=ifelse(surv_status==0,"Alive","Dead")) %>% 
    ggplot(aes(status,index_minus,fill=status))+
    geom_boxplot(alpha=0.8)+
    ggpubr::stat_compare_means(comparisons = list(c("Alive","Dead")),
                               label="p.signif")+
    theme_bw(base_size = 18)+
    theme(legend.position = "top")+
    scale_fill_manual(values= brewer.pal(n = 2, name = "Pastel2"))+
    labs(y="Pathway index",
         fill="Status")+
    x11()
  
  
  
  dat_group %>% 
    drop_na() %>% 
    ggplot(aes(age,index_minus,fill=age))+
    geom_boxplot(alpha=0.8)+
    ggpubr::stat_compare_means(comparisons = list(c("<65","≥65")),
                               label="p.signif")+
    theme_bw(base_size = 18)+
    theme(legend.position = "top")+
    scale_fill_manual(values= brewer.pal(n = 2, name = "Pastel2"))+
    labs(y="Pathway index",
         fill="Status")+
    x11()
  
  
  

  library(ggstatsplot)
  
  ggbarstats(
    data         = dat_group,
    x            = pathway_index,
    y            = stage,
    title        = "GSE30219", ## title for the plot
    legend.title = "Pathway index", 
    xlab             = "",
   # ggtheme          = hrbrthemes::theme_tinyhand(),
    ggplot.component = list(ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))),
    palette          = "Set3"
  )+x11()
  
 ggsave(pp,"pp.pdf")
}

dat_group$histo %>% table()
dat_group$ind %>% table()
library(ezcox)
dat_group <- select(dat_group,-index_minus)
mymodel_sig <- ezcox(dat_group,
                     time="surv_time",
                     status = "surv_status",
                     covariates = names(dat_group)[3:ncol(dat_group)],
                     controls = NULL,
                     global_method="wald",
                     verbose=F,
                     return_models=TRUE)


mds <- get_models(mymodel_sig)
#str(mds, max.level = 1)
x11()
show_models(mds,
            merge_models = TRUE, 
            drop_controls = F)




mymodel_mul <- ezcox(dat_group,
                     time="surv_time",
                     status = "surv_status",
                     covariates = names(dat_group)[3],
                     controls = names(dat_group)[4:ncol(dat_group)],
                     global_method="wald",
                     verbose=F,
                     return_models=TRUE)

mds <- get_models( mymodel_mul)
#str(mds, max.level = 1)
x11()
show_models(mds,merge_models = T, 
            drop_controls = F)


sig_forest
## External validation
library(IMvigor210CoreBiologies)
library(tidyverse)
# imvigor210
load("dat_input/imvirgor210nes.Rda")
IMvigor210_clin %>% 
  select(os,censOS) %>% 
  rename("surv_time"=os,
         "surv_status"=censOS)->
  imvigor_pheno.clean.os
IMvigor210_fpkm.NES %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(index_minus=mtor_up-mtor_down)->
  IMvigor210_fpkm.df
dat_imv <- Creat_dat(IMvigor210_fpkm.df,
                  imvigor_pheno.clean.os,
                  "index_minus")


dat_imv %>% 
  rownames_to_column("sample") %>% 
  inner_join(IMvigor210_clin %>% 
               rownames_to_column("sample")) %>% 
  column_to_rownames("sample") %>% 
  mutate(indgroup=ifelse(index_minus>=0.108,
                         "high","low"))->dat_imv_index



dat_imv_index$`Best Confirmed Overall Response` %>% 
  as.character() %>% 
  unique() %>% 
  na_if("NE") %>% 
  na.omit() %>% 
  combn(2,simplify = F)->mycom
library(RColorBrewer)
dat_imv_index %>% 
  filter(`Best Confirmed Overall Response`!="NE") %>%
  mutate(`Best Confirmed Overall Response`=factor(`Best Confirmed Overall Response`,
                                                  levels = c("PD","SD","PR","CR"))) %>% 
  ggplot(aes(`Best Confirmed Overall Response`,
           index_minus,fill=`Best Confirmed Overall Response`))+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.75)+
  geom_boxplot(width=0.2,
               size=0.75,outlier.colour = NA )+
  ggpubr::stat_compare_means(comparisons = mycom,
                             label = "p.signif")+
  ggpubr::stat_compare_means(label.y = 2,
                             size=5)+
  geom_hline(yintercept = 0.108,
             linetype=2,col="red",size=0.8,
             alpha=0.5)+
  theme_bw(base_size = 18)+
  scale_fill_manual(values = brewer.pal(n = 4, name = "Set3"))+
  theme(legend.position = "top")+
  labs(fill="Response",
       x="Overall Reponse",
       y="Pathway index",
       title="IMvigor210")+
   scale_y_continuous(
    breaks = c(-1.2,-1,0,0.108,1,2),
    labels = c(-1.2,-1,0,0.108,1,2))+
  x11()


library(ggstatsplot)

dat_imv_index %>% 
  filter(`Best Confirmed Overall Response`!="NE") %>% 
  mutate(`Best Confirmed Overall Response`=factor(`Best Confirmed Overall Response`,
                                                  levels = c("PD","SD","PR","CR"))) %>% 
           select(`Best Confirmed Overall Response`,
                  indgroup)->piedf
           
ggpiestats(
  data         = piedf,
  x            = indgroup,
  y            = `Best Confirmed Overall Response`,
  package      = "RColorBrewer",
  palette      = "Set2",
  title        = "IMvigor210", ## title for the plot
  legend.title = "Pathway index", ## title for the legend
  caption      = ""
  
)+
  theme_ggstatsplot()
x11()

