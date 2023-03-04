
# Expression spectrum analysis
library(tidyverse)
library(ggpubr)
load("dat_output/Mtor通路研究基因.Rdata")
load("dat_input/NSCLC-exp-clin-data.Rdata")
load("dat_input/GTex_fpkm_anno_pheno.Rdata")


NSCLC_clin <- rbind(LUAD_clin_tab,LUSC_clin_tab)
NSCLC_fpkm <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)


# Mtor gene identification and matching ------------------------------------------------------------------

str_replace(mtor_hubgene,"GRK2","ADRBK1") %>% n_distinct()

interest_gene <- sort(unique(c(mtor_hubgene,motor_gene)))
interest_gene <- sort(replace(interest_gene,interest_gene=="GRK2","ADRBK1"))

# Pathway-related genes in the literature
gene_publish <- c("PIK3R1", "PIK3CA", "PIK3C2A", 
                  "PIK3C3",
                  "AKT1", "PIK3R2", "PIK3CB", 
                  "PIK3C2B", "PIK3R4", "AKT2", 
                  "PIK3R3", "PIK3CD", "PIK3C2G", "PTEN", 
                  "AKT3", "TSC2", "MTOR","RHEB",
                  "PDK1","INPP4B","PHLPP1","STK11",
                  "TSC1")


if (F) {
  library(ggvenn)
  library(showtext)
  mylist <- list("Genes from Pubmed"=unique(c(gene_publish,  gene_new_ind)),
                 "Genes from WGCNA"=mtor_hubgene,
                 "Genes from Broad"=motor_gene)

  
  ggvenn( mylist,
         fill_color = c( "#00CBCC","#FB6501","red"),
         text_size=5)+
    theme(text=element_text(family = "serif"),
          plot.title = element_text(size=25,
                                    hjust = 0.5))+
    x11() 

 openxlsx::write.xlsx(mylist,"dat_output/基因列表.xlsx")
 unique(c(gene_publish,  gene_new_ind)) %>% 
   sort() %>% 
   as.data.frame() %>% 
   write.csv("ge.csv")
 
 
 c(unique(c(gene_publish,  gene_new_ind)),
   mtor_hubgene,
   motor_gene) %>% n_distinct()

}


##9-22 [ADD GENE] FRAP1 SPP1、EFNA2、ITGA6、ITGB8、ITGB4、IL2RB、TNC、KLF5、HRH3
gene_new <- c("MTOR" , "SPP1" , "EFNA2", "ITGA6" ,
              "ITGB8",
              "ITGB4" ,"IL2RB" ,"TNC"  , "HRH3" )
gene_new %>% 
  strsplit("、") %>% 
  unlist()->gene_new




gene_publish <- sort(gene_publish)
gene_publish %>% n_distinct()


gene.all <- sort(unique(c(interest_gene,gene_publish)))

gene.all <- sort(unique(c(interest_gene,gene_publish,gene_new)))


gene.all%in%rownames(NSCLC_fpkm) 

## Exploration in Transcriptome
NSCLC_fpkm <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)
ts <- NSCLC_fpkm[1:20]

NSCLC_fpkm %>% 
  t() %>% 
  as.data.frame() %>% 
  select(any_of(gene.all)) %>% 
  mutate(group=if_else(substr(names(NSCLC_fpkm),14,15)=="11",
                        "TCGA-normal","TCGA-cancer")) %>% 
  pivot_longer(-group,names_to = "gene",values_to = "exp")->
  NSCLC_fpkm_gene.all


# Normal lung tissue data
GTex_fpkm_anno_pheno %>% 
  mutate(across(where(is.numeric),~log2(((.x^2)-0.001)+1))) %>% 
  filter(tissue=="Lung") %>% 
  mutate(group="GTEx-normal") %>% 
  select(group,any_of(gene.all))->
  GTex_fpkm_lung

gene_new%in%names(GTex_fpkm_lung)

GTex_fpkm_lung %>% 
  pivot_longer(-group,names_to = "gene",values_to = "exp")->
  GTex_fpkm_gene.all
  

gene_exp.all <- rbind(NSCLC_fpkm_gene.all,GTex_fpkm_gene.all)

gene_exp.all.list <- group_split(gene_exp.all,gene)

 



# Batch comparison --------------------------------------------------------------------
# FUNCTION
com_ind <- function(tab){
  
  ind1 <- case_when(tab$med[1]>tab$med[2]~ paste0(tab$group[1],">",tab$group[2]),
                    tab$med[1]<tab$med[2]~ paste0(tab$group[1],"<",tab$group[2]),
                    T~paste0(tab$group[1],"=",tab$group[2])
  )
  
  ind2 <- case_when(tab$med[1]>tab$med[3]~ paste0(tab$group[1],">",tab$group[3]),
                    tab$med[1]<tab$med[3]~ paste0(tab$group[1],"<",tab$group[3]),
                    T~paste0(tab$group[1],"=",tab$group[3]))
  
  ind3 <- case_when(tab$med[2]>tab$med[3]~ paste0(tab$group[2],">",tab$group[3]),
                    tab$med[2]<tab$med[3]~ paste0(tab$group[2],"<",tab$group[3]),
                    T~paste0(tab$group[2],"=",tab$group[3]))
  
  c(ind1,ind2,ind3)
  
}

compare <- function(gene_dat){
  
#Produce comparison results
gene_dat%>% 
  group_by(group) %>% 
  summarise(med=median(exp)) %>% 
  arrange(desc(group)) %>% 
  com_ind()->
  compare_ind

#Generate comparative statistical results
compare_means(exp~group,gene_dat) %>% 
 mutate(gene=gene_dat$gene[1],compare=compare_ind) %>% 
  select(gene,everything()) %>% 
  select(-2)->
  res

res
}


gene_dat<- gene_exp.all.list[[1]]




gene_exp.all.list %>% 
  map_df(compare)->
  res_all


setwd("C:/Users/Cloudy/Desktop/论文写作/文章写作预分析/dat_output")
save(res_all,file = "兴趣基因批量肿瘤与正常样本表达对比.Rdata")
load("dat_output/兴趣基因批量肿瘤与正常样本表达对比.Rdata")

#Gene filtration--------------------------------------------------------------------
#We hope to include those genes with the following expression patterns
##1. Tumor ＜ paracancerous ≤ normal
##2. Tumor ＜ paracancerous ≥ normal (not considered temporarily)
##3. Tumor>para-cancer ≥ normal
##4. Tumor ＞ para-cancer ≤ normal (not considered temporarily)
##View actual groups

res_all$compare %>% 
  unique() %>% 
  sort()

## Build Group



load("dat_output/兴趣基因批量肿瘤与正常样本表达对比.Rdata")

## Build function
check_gene <- function(tab){
  
  #Compare pairs
  mycom <- tab %>% pull(compare)
  # GENE
  gene <- tab %>% pull(gene) %>% unique()
  # Compare the judgment of existence
  com1 <- all(c("TCGA-normal>TCGA-cancer",
                "TCGA-cancer<GTEx-normal",
                "TCGA-normal<GTEx-normal")%in%mycom)
  
  com2<- all(c("TCGA-normal<TCGA-cancer",
               "TCGA-normal>GTEx-normal",
               "TCGA-cancer>GTEx-normal")%in%mycom)
  

  check_res <- ifelse(any(c(com1,com2)),gene,NA) 
  
}

res_all %>% 
 group_split(gene) %>% 
 map_chr(check_gene) %>% 
 na.omit() %>% 
 as.character()->
 gene_compare.satis


## The effect of gene on survival is not considered, only the effect of expression is considered

res_all %>% 
  filter(compare=="TCGA-normal<TCGA-cancer"&p.signif=="ns"|
           compare=="TCGA-normal>TCGA-cancer"&p.signif=="ns"| 
           compare=="TCGA-normal>GTEx-normal"&p.signif=="ns"| 
           compare=="TCGA-normal<GTEx-normal"&p.signif=="ns"|  
           compare=="TCGA-cancer<GTEx-normal"&p.signif=="ns"|  
           compare=="TCGA-cancer>GTEx-normal"&p.signif=="ns") %>% 
  pull(gene)->
  gene.exp.ns #No difference in the expression of 66 genes



##Genes meeting the required expression trend
res_all %>% 
  filter(gene%in%gene_compare.satis) %>%  #Meet expression mode
  filter(!gene%in%gene.exp.ns) %>%   #Statistically significant
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
  gene_exp.sig

gene_exp.sig %>% 
  group_split(gene) %>% 
  map_df(~mutate(.,group_com=paste(.$compare,collapse = "/"))) %>% 
  pull(group_com) %>% 
  table()



gene_exp.sig$gene %>% n_distinct()

gene_new%in%gene_exp.sig

gene_exp.sig %>%
  pull(gene) %>%
  n_distinct()#69 genes meet the group screening requirements

gene_new[gene_new%in%gene_exp.sig$gene]->gene_expnew

intersect(c(gene_newos,gene_newdfs),gene_expnew)

# New inclusion genes
gene_new_ind <- setdiff(gene_new,"KLF5")




os.res %>% 
  filter(gene!="KLF5") %>% 
  pull(gene)->gene1



"MTOR"%in%rownames(NSCLC_fpkm)


gene_exp.sig %>% 
  openxlsx::write.xlsx("满足表达同相性的通路基因.xlsx")

## Consider the impact of genes on survival
load("dat_output/通路基因生存分析结果.Rdata")
list("OS of pathway related genes"=os.res,
     "DFS of pathway related genes"=dfs.res) %>% 
  openxlsx::write.xlsx("dat_output/通路相关基因有生存的影响结果.xlsx")

c(gene_danger.all[gene_danger.all%in%gene_compare.satis],
gene_protect.all[gene_protect.all%in%gene_compare.satis])->
  gene_fin


## Filter out the genes that have no difference in the expression of tumor and tumor
res_all %>% 
  filter(gene%in%gene_fin) %>%
  mutate(group=ifelse(gene%in%gene_protect.all,"gene.protect","gene.danger")) %>% 
  filter(compare=="TCGA-normal<TCGA-cancer"&p.signif=="ns"|
         compare=="TCGA-normal>TCGA-cancer"&p.signif=="ns"| 
         compare=="TCGA-normal>GTEx-normal"&p.signif=="ns"| 
         compare=="TCGA-normal<GTEx-normal"&p.signif=="ns"|  
         compare=="TCGA-cancer<GTEx-normal"&p.signif=="ns"|  
         compare=="TCGA-cancer>GTEx-normal"&p.signif=="ns"|   
         compare=="TCGA-normal<TCGA-cancer"&group=="gene.protect"|
         compare=="TCGA-normal>TCGA-cancer"&group=="gene.danger") %>% 
  pull(gene)->
  gene.ns






res_all %>% 
  filter(gene%in%gene_fin) %>%
  filter(!gene%in%gene.ns) %>% 
  mutate(group=ifelse(gene%in%gene_protect.all,"gene.protect","gene.danger")) %>% 
  pull(gene) %>% 
  unique()-> 
  gene.important





## Gene mapping
gene.important%in%gene_protect.both %>% sum()
#2
gene.important%in%gene_danger.both%>% sum()
#1
gene.important%in%os_res_danger$gene%>% sum()
#11
gene.important%in%os_res_protect$gene%>% sum()
#9
gene.important%in%dfs_res_danger$gene%>% sum()
#2
gene.important%in%dfs_res_protect$gene%>% sum()
#4



res_all %>% 
  filter(gene%in%gene.important) %>% 
  inner_join(os.res %>% 
               mutate(index="OS") %>% 
               select(gene,index,HR,`Cox p`, `Log-rank p`) %>% 
               bind_rows(dfs.res %>% 
                       mutate(index="DFS") %>% 
                      select(gene,index,HR,`Cox p`, `Log-rank p`))) %>% 
  arrange(index)->
  gene.important.tab
names(os.res)
# Preserve important genes
save(gene.important,
     gene_protect.all,
     gene_protect.both,
     gene.important.tab,
     gene_danger.all,gene_danger.both,
     file = "dat_output/通路基因生存分析及表达谱分析结果.Rdata")

