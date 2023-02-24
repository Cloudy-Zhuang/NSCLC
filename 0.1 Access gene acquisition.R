setwd("C:/Users/Cloudy/Desktop/论文写作/文章写作预分析/dat_input")
library(tidyverse)
library(data.table)

#Data download from: https://www.kegg.jp/kegg/rest/keggapi.html
#Consult：https://zhuanlan.zhihu.com/p/434383719
hsa_pathway <- fread("hsa.txt",data.table = F,header = F)
hsa_gene <- fread("hsa_gene.txt",data.table = F,header = F)
path_pathway <- fread("path_pathway.txt",data.table = F,header = F)

hsa_pathway %>% 
  mutate(V2=str_remove_all( V2," - Homo sapiens \\(human\\)")) %>% 
  mutate_all(trimws) %>% 
  set_names("path","pathway.name") %>% 
  inner_join(path_pathway %>% 
               set_names("path","ID")) %>%  
  select(-path)->
  hsa_pathway_clean


hsa_gene %>% 
  separate(V2,into = c("gene","fullname"),sep = ";") %>%
  mutate(genecopy=gene) %>% 
  separate_rows(gene,sep = ",") %>% 
  set_names("ID","gene","fullname","genecopy") %>% 
  mutate_all(trimws)->
  hsa_gene_clean

pathwaydata <- inner_join(hsa_pathway_clean,
                          hsa_gene_clean)


## Extract PI3K path information
pathwaydata %>% 
  filter(str_detect(pathway.name,"PI3K"))->PI3K

#save data
setwd("C:/Users/Cloudy/Desktop/论文写作/文章写作预分析/dat_output")

save(PI3K,file="pi3k_inf.Rdata")

