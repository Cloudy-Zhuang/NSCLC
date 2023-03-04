library(tidyverse)
library(ggpubr)
load("dat_output/通路基因生存分析及表达谱分析结果.Rdata")
load("dat_input/NSCLC-exp-clin-data.Rdata")

NSCLC_fpkm <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)
NSCLC_clin <- rbind(LUAD_clin_tab,LUSC_clin_tab)



#Tumor - paracancerous - normal lung tissue-------------------------------------------------------------
#res_ All result sorting


res_all %>% 
  mutate(compare=str_replace_all(compare,
                "TCGA-normal<TCGA-cancer",   #Adjust the order
                  "TCGA-cancer>TCGA-normal")) %>% 
  mutate(compare=str_replace_all(compare,
                               "TCGA-normal>TCGA-cancer",
                               "TCGA-cancer<TCGA-normal")) %>% 
  group_split(gene) %>% 
  map_df(~arrange(.x,compare)) %>% 
  group_split(gene) %>% 
  map_df(~mutate(.x,compare=c(.x$compare[c(3,2,1)]))) %>% 
  group_split(gene) %>% 
  map_df(~data.frame(Gene=.x$gene[1],
                     Group1=paste0(.x$compare[1]," (","p.signif: ",.x$p.signif[1],")"),
                     Group2=paste0(.x$compare[2]," (","p.signif: ",.x$p.signif[2],")"),
                     Group3=paste0(.x$compare[3]," (","p.signif: ",.x$p.signif[3],")")
                     ))->
  res_all.adjust  #Output an experimental table




res_all %>% 
  mutate(compare=str_replace_all(compare,
                                 "TCGA-normal<TCGA-cancer",
                                 "TCGA-cancer>TCGA-normal")) %>% 
  mutate(compare=str_replace_all(compare,
                                 "TCGA-normal>TCGA-cancer",
                                 "TCGA-cancer<TCGA-normal")) %>% 
  group_split(gene) %>% 
  map_df(~arrange(.x,compare)) %>% 
  group_split(gene) %>% 
  map_df(~mutate(.x,compare=c(.x$compare[c(3,2,1)])))->
  res_sort

##  Grouping
res_sort %>% 
  pull(compare) %>% 
  unique()->
  compare_group



res_sort %>%
  filter(p<0.05) %>% 
  group_split(compare) %>% 
  map_chr(~.x$compare[1])->
  groupname

res_sort %>%
  filter(p<0.05) %>% 
  select(gene,compare) %>% 
  group_split(compare) %>% 
  map(~select(.x,-compare)) %>% 
  map(unlist) %>% 
  map(`length<-`,max(map_int(.,length))) %>% 
  bind_cols() %>% 
  set_names(groupname)->
  gene_compare
  


res_sort %>%
  filter(p<0.05) %>% 
  select(gene,compare) %>% 
  group_split(compare) %>% 
  map(~select(.x,-compare)) %>% 
  map(unlist) %>% 
  map_int(length) %>%
  data.frame() %>% 
  mutate(group=groupname) %>% 
  set_names("Count","group")->
  exp_com




###Save file

setwd("C:/Users/Cloudy/Desktop/论文写作/文章写作预分析/table")
res_sort %>% 
  select(-p.adj) %>%
  rename("Comparison of median expression levels"=
           compare) %>% 
  set_names(str_to_title(names(.))) %>% 
  openxlsx::write.xlsx("基因在不同组织中的表达对比.xlsx")
gene_compare %>% 
  openxlsx::write.xlsx("不同组织中对比有差异的基因.xlsx")





# Introduction of mutation data ------------------------------------------------------------------
# Ref: https://zhuanlan.zhihu.com/p/419960451

load("dat_output/TCGA突变数据.Rdata")

 
# NSCLC interest mutation data processing

NSCLC_sen.mut %>% 
  filter(substr(sample,1,12)%in%substr(names(NSCLC_fpkm),1,12)) %>% 
  mutate(across(-sample,~ifelse(is.na(.x),"wild","mutation")))->
  NSCLC_interest_mut





NSCLC_fpkm %>% 
  t() %>% 
  as.data.frame() %>% 
  select(gene.important) %>% 
  mutate(group=ifelse(substr(rownames(.),14,15)=="11","normal","cancer"))->
  roc.pre.dat
# save(roc.pre.dat,file="dat_output/roc准备数据.Rdata")

load("dat_output/mtor通路基因清洁版本.Rdata")


# Extract variables with missing values less than two-thirds except OS
NSCLC_clin %>%
  filter(rownames(.)%in%substr(names(NSCLC_fpkm),1,12)) %>% 
  map_lgl(~sum(is.na(.x))<597)->logind




interest_ind <- names(NSCLC_clin)[logind][c(1,3:9,11)]

## Build calculated functions
if (F) {
  Compare_means<- function(gene, dat) {
    # Calculate the median value
    median_res <- sort(tapply(dat[,gene],dat$group,median))
    med_compare <-case_when(median_res[1]>median_res[2]~paste0(names(median_res[1])," > ",names(median_res[2])),
                             median_res[1]<median_res[2]~paste0(names(median_res[1])," < ",names(median_res[2])),
                             T~paste0(names(median_res[1]),"=",names(median_res[2]))) 
    
    # Homogeneity test of variance
    var_check <- var.test(dat[, gene] ~ group, data = dat)
    # If the variance is 0 or cannot be calculated as the missing value, the data will not be able to perform the normal distribution test, so the missing value is directly set and returned
    if (var_check$estimate == 0 | is.na(var_check$estimate)) {
      
      return(data.frame(gene = gene, statistic = NA, pvalue = NA, method = NA),
             median.compare=med_compare)
    }
    # Normality test
    shap <- tapply(dat[, gene], dat$group, shapiro.test)
    
    # If any one of the two groups of data conforms to the normal distribution and the comparison variance of the two groups of data is the same, perform the t-test, otherwise execute the wilcox.test
    # Determine whether to perform t test or rank sum test and return statistics and P value
    if (all(c(shap[[1]]$p.value, shap[[2]]$p.value) > 0.05 && var_check$p.value > 0.05)) {
      test <- t.test(dat[, gene] ~ group, data = dat)
      res <- data.frame(gene = gene, statistic = test$statistic, 
                        pvalue = test$p.value, 
                        method = "t.test",
                        median.compare=med_compare )
    } else {
      test <- wilcox.test(dat[, gene] ~ group, data = dat, exact = F)
      res <- data.frame(gene = gene, statistic = test$statistic, 
                        pvalue = test$p.value, method = "wilcox.test",
                        median.compare=med_compare )
    }
    
    return(res)
  }
  
}


##Pathway gene expression matrix

NSCLC_fpkm %>% 
  t() %>% 
  as.data.frame() %>% 
  select(pathway.gene.all) %>% 
  mutate(sample=substr(rownames(.),1,12)) %>% 
  select(sample,everything()) ->
  mtorgene.exp

NSCLC_clin$Smoking_exposure





## Build batch operation function

if (F) {
multicom <-  function(ind){
  print(ind)
##Combination of clinical factors and expression matrix
    NSCLC_clin %>% 
      mutate(sample=rownames(.)) %>% 
      select(sample,ind) %>% 
      inner_join(mtorgene.exp) %>% 
      distinct(sample,.keep_all = T) %>% 
      select(-sample) %>% 
      select(ind,everything()) %>% 
      na_if("") %>% 
      drop_na() %>% 
      rename("group"=ind)->
      expdat
    
    map_df(names(expdat)[-1],
           Compare_means,expdat) %>% 
      mutate(index=ind) ->res
    res
    
}

}
# Batch operation
clin_compare_res<- lapply(interest_ind,multicom)
#save(clin_compare_res ,file = "dat_output/通路基因在不同临床信息分组中的表达情况.Rdata")


# Mutation --------------------------------------------------------------------
## Mutations are targeted at tumor samples, so they should not be combined with clinical matrix
ind <- names(NSCLC_interest_mut)[2]

# Build calculated functions
if (F) {
  mutcom <- function(ind){
    
    mtorgene.exp %>% 
      filter(substr(rownames(.),14,15)!="11") %>%  #Filter non-tumor samples
      mutate(sample=substr(rownames(.),1,15)) %>% 
      inner_join(NSCLC_interest_mut %>% 
                   select(sample,ind)) %>% 
      select(-sample) %>%
      rename("group"=ind) %>% 
      select(group,everything())->
      mut_dat
    
    map_df(names(mut_dat)[-1],
           Compare_means,mut_dat) %>% 
      mutate(index=ind) ->res
    res
  }
}

mut_compare_res<- lapply(names(NSCLC_interest_mut)[-1],mutcom)

#save(mut_compare_res,file = "dat_output/通路基因在不同兴趣突变分组中的表达.Rdata")
# Result summary and output -----------------------------------------------------------------

rm(list = ls())

load("dat_output/通路基因在不同临床信息分组中的表达情况.Rdata")
load("dat_output/通路基因在不同兴趣突变分组中的表达.Rdata")

tab.res_attain <- function(tab){

 tab %>% 
  filter(pvalue<0.05)->
  sig.tab

 sig.tab$median.compare %>% 
  table() %>% # Count the number of genes in different groups
  as.data.frame() %>%
  set_names("group","Freq") %>% 
  mutate(Count=sum(.$Freq)) %>%  #Calculate the total number of genes
  mutate(index=unique(sig.tab$index)) %>% 
  select(index,everything()) %>% 
  arrange(index)->
  res.tab
 
 res.tab
}

clin.count.tab <- map_df(clin_compare_res,tab.res_attain)
mut.count.tab <- map_df(mut_compare_res,tab.res_attain)
clin.count.tab %>% 
  DT::datatable()
mut.count.tab %>% 
  DT::datatable()
## Gene results and table output
setwd("C:/Users/Cloudy/Desktop/论文写作/文章写作预分析/table")


clin_compare_res %>% 
  set_names(clin_compare_res %>% 
              map_chr(~pull(.x,index) %>% 
                        unique)) %>% 
  map(~filter(.x,pvalue<0.05) %>% 
        select(-2) %>% 
        rename("Comparison of median expression values"=4)) %>% 
  openxlsx::write.xlsx("不同临床性状中表达有差异的基因.xlsx",
                       over.wirte=T)


mut_compare_res %>% 
  set_names(mut_compare_res %>% 
              map_chr(~pull(.x,index) %>% 
                        unique)) %>% 
  map(~filter(.x,pvalue<0.05) %>% 
        select(-2) %>% 
        rename("Comparison of median expression values"=4)) %>% 
  openxlsx::write.xlsx("不同突变表型中表达有差异的基因.xlsx")
  



exp.comres <- openxlsx::read.xlsx("基因在不同组织中的表达对比.xlsx")
exp.interest <- openxlsx::read.xlsx("满足表达同相性的通路基因.xlsx")


exp.interest %>% 
  group_split(gene) %>% 
  map_df(~data.frame(gene=unique(.x$gene),
                     group=paste0(.x$compare[1],"/",
                                  .x$compare[2],"/",
                                  .x$compare[3] )))->
  exp.combine


exp.combine %>% 
  pull(group) %>% 
  table() %>% 
  as.data.frame() %>% 
  set_names('group',"Freq") %>% 
  mutate(Count=sum(.$Freq)) %>% 
  mutate(index="Expression") %>% 
  relocate("index")->
  expcom.tab


exp.comres %>% 
  filter(P<0.05) %>% 
  pull(8) %>% 
  table() %>% 
  as.data.frame() %>% 
  set_names("group","Freq") %>% 
  mutate(Count=exp.comres %>% 
           filter(P<0.05) %>% 
           pull(Gene) %>% 
           n_distinct()) %>% 
  mutate(index="Expression") %>% 
  relocate(index)->
  exp.count.tab
exp.count.tab %>% 
  DT::datatable()

tab.all <- rbind(exp.count.tab,
                 expcom.tab,
                 clin.count.tab,
                 mut.count.tab
                 )

tab.all %>% 
  DT::datatable()

tab.all %>% 
  openxlsx::write.xlsx("通路基因与多性状表达汇总.xlsx")
  

# Find a gene with no difference in expression
# exp.comres %>%
#   filter(P>0.05) %>%
#   pull(Gene) %>%
#   unique()->gene.ns
# 
# exp.comres %>%
#   filter(P<0.05) %>%
#   pull(Gene) %>%
#   unique()->gene.sig
# 
# "CDKN1B"
# setdiff(gene.ns,gene.sig)






