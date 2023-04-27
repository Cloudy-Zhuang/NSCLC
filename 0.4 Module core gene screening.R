#Module core gene screening
library(tidyverse)


#Read in data
load("dat_output/Preliminary results of core genes in pi3k_akt_mtor_turquoise module after module merging.Rdata")


pathway_gene <- clusterProfiler::read.gmt("dat_input/HALLMARK_PI3K_AKT_MTOR_SIGNALING.v7.5.1.gmt") %>% 
                pull(gene)

Find_hubgene <- function(pathway_res,cancer_res){
  
pathway_res %>% 
    filter(abs(MM)>0.8&abs(GS)>0.2&MMP<0.05&GSP<0.05) %>% 
    pull(moduleGenes)->gene_path
    
cancer_res %>% 
    filter(abs(MM)>0.8&abs(GS)>0.2&MMP<0.05&GSP<0.05) %>% 
    pull(moduleGenes)->gene_dis
  
  gene_res <- intersect(gene_path,gene_dis)
  gene_res
}



hubgene <- Find_hubgene(mydata_mtor,mydata_mtor_cancer)
intersect(hubgene,pathway_gene) 

mtor_hubgene <- hubgene
motor_gene <- pathway_gene

save(mtor_hubgene,
     motor_gene,
     file="dat_output/The PI3K AKT Mtor pathway genes for study.Rdata")


