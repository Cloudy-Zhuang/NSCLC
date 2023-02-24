library(tidyverse)
library(clusterProfiler)
setwd("C:/Users/Cloudy/Desktop/论文写作/文章写作预分析")
load("dat_input/NSCLC_disgene.Rdata")
load("dat_input/TCGA_NSCLC_DEGs.Rdata")
load("dat_output/pi3k_inf.Rdata")



PI3K %>%
  mutate(genename.con=str_split_fixed(genecopy,",",2)[,1])->
  PI3K
  
TCGA_NSCLC_diffgene %>% 
  filter(padj<0.05&abs(log2FoldChange)>=1) %>% 
  filter(gene_id%in%PI3K$gene) %>% 
  inner_join(PI3K %>% 
               select(gene,fullname) %>% 
               rename("gene_id"=gene))-> pi3k_genediff
  

pi3k.akt <- read.gmt("dat_input/WP_PI3KAKT_SIGNALING_PATHWAY.v7.5.1.gmt")

TCGA_NSCLC_diffgene %>% 
  filter(padj<0.05&abs(log2FoldChange)>=1) %>% 
  filter(gene_id%in%pi3k.akt$gene) %>% 
  inner_join(PI3K %>% 
               select(gene,fullname) %>% 
               rename("gene_id"=gene))-> pi3k_genediff_office


pi3k_genediff$gene_id %>% 
  setdiff(pi3k_genediff_office$gene_id)

PI3K %>% 
  filter(gene%in%c("ERBB4" ,"FIGF" ))



# save data
save(pi3k_genediff,file="PI3K通路中的差异基因.Rdata")



# mtor path matching ----------------------------------------------------------------

mtor <- read.gmt("dat_input/HALLMARK_PI3K_AKT_MTOR_SIGNALING.v7.5.1.gmt")

TCGA_NSCLC_diffgene %>% 
  filter(padj<0.05&abs(log2FoldChange)>=1) %>% 
  filter(gene_id%in%mtor$gene ) %>% 
  inner_join(PI3K %>% 
               select(gene,fullname) %>% 
               rename("gene_id"=gene))-> mtor_genediff
setwd("C:/Users/Cloudy/Desktop/论文写作/文章写作预分析/dat_output")

save(mtor_genediff,file="mtor通路差异基因.Rdata")

