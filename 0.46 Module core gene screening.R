#Module core gene screening
library(tidyverse)


#Read in data
load("dat_output/合并模块后pi3k_akt_mtor_turquoise模块核心基因初步结果.Rdata")
load("dat_output/pi3k_akt.tan模块核心基因初步结果.Rdata")

pathway_gene <- clusterProfiler::read.gmt("dat_input/HALLMARK_PI3K_AKT_MTOR_SIGNALING.v7.5.1.gmt") %>% 
                pull(gene)

pathway_gene <- clusterProfiler::read.gmt("dat_input/WP_PI3KAKT_SIGNALING_PATHWAY.v7.5.1.gmt") %>% 
  pull(gene)

#ref: https://mp.weixin.qq.com/s/tQ0wLO74I80-DmgGy2mC_A

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
hubgene <- Find_hubgene(merge_pi3k_akt,merge_pi3k.akt_cancer)
intersect(hubgene,pathway_gene) 

mtor_hubgene <- hubgene
motor_gene <- pathway_gene

save(mtor_hubgene,motor_gene,file="Mtor通路研究基因.Rdata")


