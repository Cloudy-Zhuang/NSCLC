# Gene expression profile analysis: combined with GTex
setwd("E:/突击文章/7.GTex 数据")
library(tidyverse)
library(data.table)

# GTex data cleaning ---------------------------------------------------------------

dat <- fread("gtex_RSEM_gene_fpkm.gz",data.table = F)
anno <- fread("probeMap_gencode.v23.annotation.gene.probemap",data.table = F)
pheno <- fread("GTEX_phenotype.gz",data.table = F)

dat %>% 
  rename("id"=sample) %>% 
  inner_join(select(anno,id,gene)) %>% 
  select(-id) %>% 
  select(gene,everything()) %>% 
  mutate(ind=rowSums(across(where(is.numeric)))) %>% 
  arrange(desc(ind)) %>% 
  distinct(gene,.keep_all = T) %>% 
  select(-ind) %>% 
  column_to_rownames("gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") ->GTex_fpkm_anno
  

GTex_fpkm_anno %>% 
  inner_join(select(pheno,Sample,`_primary_site`)) %>% 
  select(Sample,`_primary_site`,everything()) %>% 
  rename("tissue"=`_primary_site`)->GTex_fpkm_anno_pheno
setwd("E:/突击文章/dat_out")
#save(GTex_fpkm_anno_pheno,file = "GTex_fpkm_anno_pheno.Rdata")




# Generate a Getx matrix of interest genes -----------------------------------------------------------
rm(list = ls())
setwd("E:/突击文章/dat_out")
load("GTex_fpkm_anno_pheno.Rdata")
load("Modelgenes.Rdata")
GTex_fpkm_anno_pheno$tissue %>% table()
GTex_fpkm_anno_pheno %>% 
  filter(tissue%in%"Lung") %>% 
  select(tissue,strategy_sec_genes)->GTex_modelgenedata
save(GTex_modelgenedata,file = "GTex_modelgenedata.Rdata")
# Interest genes and TCGA-Fpkm data were introduced -----------------------------------------------------
load("Modelgenes.Rdata")
load("NSCLC-COVID-drug单因素cox分析结果.Rdata")
load("NSCLC-exp-clin-data.Rdata")
load("GTex_modelgenedata.Rdata")

# Representation matrix merging ------------------------------------------------------------------

NSCLC_fpkmanno <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)
test <- NSCLC_fpkmanno[1:20]

rownames(ts) %>% 
  substr(14,15) %>% 
  table()


NSCLC_fpkmanno %>% 
  filter(rownames(.)%in%cox_res$gene) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(group=if_else(substr(rownames(.),14,15)=="11",
                       "TCGA-Normal","TCGA-NSCLC")) %>% 
  dplyr::select(group,cox_res$gene,everything()) ->NSCLC_coxgene 


GTex_modelgenedata %>% 
  dplyr::select(cox_res$gene) %>% 
  mutate(group="GTEX-Normal") %>% 
  bind_rows(NSCLC_coxgene)->modelgene_expression
  


# draw ----------------------------------------------------------------------

dat <- modelgene_expression
mycom <- combn(unique(dat$group),2,simplify = F)
mylist <- list()
i=6
dat <- dplyr::select(dat,sort(names(dat)[-ncol(dat)]),group)

 #myplot <- list()
for( i in 1:(ncol(dat)-1)) {
  
  ymax=max(dat[[i]])*1.8
  ymin=min(dat[[i]])
  step = 0.08
# myplot[[i]]
myplot<- ggplot(aes_string("group",names(dat[i]),fill="group"),
                      data = dat)+
    geom_boxplot(outlier.colour = NA,size = 0.4)+
    scale_fill_manual(values = c("#00ba38","#619cff","#f8766d"))+
    theme_bw(base_size = 15)+
    labs(y=paste0(names(dat[i])),x="",fill="group")+
    theme(legend.position = "none",
          #panel.spacing.x = unit(0.2,"cm"),
          #panel.spacing.y = unit(0.1, "cm"),
          strip.text.x = element_text(size=5,color="black"),
          #strip.background.x = element_blank(),
          axis.text = element_text(color="black"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text = element_text(color="black",size=11),
          #legend.title=element_blank(),
          legend.spacing.x=unit(0.1,'cm'),
          legend.key=element_blank(),
          legend.key.width=unit(0.5,'cm'),
          legend.key.height=unit(0.5,'cm')
          #plot.margin=unit(c(0.3,0.3,0.3,0.3),units=,"cm")
    )+
    ggpubr::stat_compare_means(comparisons = mycom,
                               label = "p.signif")+
  ylim(ymin-0.1, ymax+0.1*(ymax-ymin)+step*(ymax-ymin)*length(mycom))

ggsave(myplot,filename = paste0(names(dat)[i],".pdf"),height = 4,width = 4)

}

modelgene_expression %>% 
  filter(!str_detect(group,"GTEX")) ->ts




wilcox.test(ts$SLCO1B1~ts$group)


modelgene_expression$SLCO1B1 %>% summary()
modelgene_expression$SLCO1B1 %>% table()
dir.create("模型基因表达量箱图")

table(modelgene_expression$group)


## Jigsaw puzzle
library(cowplot)
x11()
a <- plot_grid(plotlist = myplot,ncol = 3,labels = letters[1:7])

getwd()


# GEO ---------------------------------------------------------------------
load("GEO纳入研究的表达矩阵.Rdata")

GSE147507_include_exp[cox_res$gene,] %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate_all(~log2(.x+1)) %>% 
  mutate(group=ifelse(grepl("SARS-CoV|COVID",rownames(.)),"COVID-19","normal"))->
  GSE147507_modelgene

GSE157103_include_exp %>%
  dplyr::rename("gene"=1) %>% 
  filter(gene%in%cox_res$gene) %>% 
  column_to_rownames("gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate_all(~log2(.x+1)) %>% 
  mutate(group=ifelse(grepl("^COVID",rownames(.)),"COVID-19","normal"))->
  GSE157103_modelgene


GSE166190_include_exp[cox_res$gene,] %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate_all(~log2(.x+1)) %>% 
  mutate(group=ifelse(grepl("Pos",rownames(.)),"COVID-19","normal"))->
  GSE166190_modelgene



dat=GSE166190_modelgene
for( i in 1:(ncol(dat)-1)) {
  
  ymax=max(dat[[i]])*1.8
  ymin=min(dat[[i]])
  step = 0.08
  
  myplot[[i]] <- ggplot(aes_string("group",names(dat[i]),fill="group"),data = dat)+
    geom_boxplot(outlier.colour = NA,size = 0.4)+
    scale_fill_manual(values = c("#00ba38","#619cff","#f8766d"))+
    theme_bw(base_size = 15)+
    labs(y=paste0(names(dat[i])),x="",fill="group")+
    theme(
      #panel.spacing.x = unit(0.2,"cm"),
      #panel.spacing.y = unit(0.1, "cm"),
      strip.text.x = element_text(size=5,color="black"),
      #strip.background.x = element_blank(),
      axis.text = element_text(color="black"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.text = element_text(color="black",size=11),
      #legend.title=element_blank(),
      legend.spacing.x=unit(0.1,'cm'),
      legend.key=element_blank(),
      legend.key.width=unit(0.5,'cm'),
      legend.key.height=unit(0.5,'cm')
      #plot.margin=unit(c(0.3,0.3,0.3,0.3),units=,"cm")
    )+
    ggpubr::stat_compare_means(label.x = 0.5,label.y = 0.5*ymax)+
    ylim(ymin-0.1, ymax+0.1*(ymax-ymin)+step*(ymax-ymin)*length(mycom))
  
  #ggsave(myplot,filename = paste0(names(dat)[i],".pdf"),height = 4,width = 5)
  
}

# Fix weird genes -----------------------------------------------------------------
library(ggpubr)
gene="GCGR"
gene_boxplot <- function(gene){
  
  ymax=max(dat[gene])*1.8
  ymin=min(dat[gene])
  step = 0.08
    
  
modelgene_expression %>% 
  filter(!str_detect(group,"GTEX")) %>% 
  ggplot(aes_string("group",gene,fill="group"))+
  geom_boxplot()+
  scale_fill_manual(values = c("#619cff","#f8766d"))+
  theme_bw(base_size = 15)+
  labs(y=gene,x="",fill="group")+
  ylim(0,0.6)+
  theme(
    #panel.spacing.x = unit(0.2,"cm"),
    #panel.spacing.y = unit(0.1, "cm"),
    strip.text.x = element_text(size=5,color="black"),
    #strip.background.x = element_blank(),
    axis.text = element_text(color="black"),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.text = element_text(color="black",size=11),
    #legend.title=element_blank(),
    legend.spacing.x=unit(0.1,'cm'),
    legend.key=element_blank(),
    legend.key.width=unit(0.5,'cm'),
    legend.key.height=unit(0.5,'cm'),
    #legend.position = "none",
    #plot.margin=unit(c(0.3,0.3,0.3,0.3),units=,"cm")
  )+stat_compare_means(comparisons = list(c("TCGA-NSCLC" ,"TCGA-Normal")),
                       label.y = 0.5,
                        label = "p.signif")#->gene_plot
  # ylim(ymin-0.9, ymax+0.1*(ymax-ymin)+
  #        step*(ymax-ymin)*1)

#ggsave(gene_plot,filename = paste0("箱图修",gene,".pdf"),height = 4,width = 4)  

}
x11()
gene_boxplot("GCGR")
gene_boxplot("IL2")
gene_boxplot("SLCO1B1")
three_gene <- lapply(list("GCGR","IL2","SLCO1B1"),gene_boxplot)

x11()
plot_grid(plotlist = three_gene,labels = letters[1:3])

# save data ---------------------------------------------------------------
setwd("C:/Users/Cloudy/Desktop/突击文章/dat_out")
save(modelgene_expression,file = "模型基因表达矩阵.Rdata")
