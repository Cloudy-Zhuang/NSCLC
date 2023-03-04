library(tidyverse)
library(openxlsx)
setwd("C:/Users/Cloudy/Desktop/论文写作/文章写作预分析/dat_output")
load("pi3k-akt-mtor通路基因生存分析结果.Rdata")


os.res %>% 
  filter(HR>1)->os_res_danger
os.res %>% 
  filter(HR<1)->os_res_protect

dfs.res %>% 
  filter(HR>1)->dfs_res_danger
dfs.res %>% 
  filter(HR<1)->dfs_res_protect

intersect(os_res_danger$gene,dfs_res_danger$gene)->gene_danger.both

intersect(os_res_protect$gene,dfs_res_protect$gene)->gene_protect.both


c(os_res_danger$gene,dfs_res_danger$gene)->gene_danger.all

c(os_res_protect$gene,dfs_res_protect$gene)->gene_protect.all

save(list=ls(pattern = "os|dfs"),file="dat_output/通路基因生存分析结果.Rdata")


