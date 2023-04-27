library(tidyverse)
library(clusterProfiler)
library(openxlsx)

## Read data
load("dat_input/NSCLC-exp-clin-data.Rdata")
gene.important <- read.xlsx("dat_input/gene.important.xlsx") %>% 
  pull(1)

pi3k_mtor <- read.gmt("dat_input/HALLMARK_PI3K_AKT_MTOR_SIGNALING.v7.5.1.gmt")
hallmark <- read.gmt("dat_input/h.all.v7.5.1.symbols.gmt")
NSCLC_fpkm <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)

NSCLC_fpkm %>% 
  t() %>% 
  as.data.frame()->NSCLC_fpkm

gene_input <- gene.important[1]
dat <- select(NSCLC_fpkm,gene_input,1:50)
# Batch dependency partition function ---------------------------------------------------------------
cor.bath <- function(gene_input,dat){
  print(gene_input)
  
  gene_exp <- dat[[gene_input]]
  
  dd <- do.call(rbind,lapply(setdiff(names(dat),
                                     gene_input),
                            function(x){
                              
  cor_res <- cor.test(as.numeric(dat[,x]),
                  gene_exp,
                  method ="spearman",
                  exact=FALSE)
  
  data.frame(gene1=gene_input,gene2=x,cor= cor_res$estimate,p.value=cor_res$p.value)}
  ))
  
 
}
ts <- cor.bath(gene_input,dat)

library(future)
library(future.apply)
plan(multisession)
res.list <- future_lapply(gene.important,cor.bath,NSCLC_fpkm)


save(res.list,file = "dat_output/gene.imp cor bath.Rdata")

load("dat_output/gene.imp cor bath.Rdata")

# GSEA analysis
 ------------------------------------------------------------------



Count_gsea <- function(tab,gmt){

gene <- tab$gene2

## geneList trilogy
## 1.Get gene logFC
geneList <- tab$cor
## 2.name
names(geneList) = gene
## 3.Sorting is important
geneList = sort(geneList, decreasing = TRUE)

# You need Internet. It's crowded, but it's fast
y <- GSEA(geneList,
          TERM2GENE =gmt,
          nPermSimple = 10000)

res <- y@result %>% 
  mutate(gene=unique(tab$gene1)) %>% 
  dplyr::select(gene,everything())

res
}

gsea_res <- map_df(res.list,
                   Count_gsea,
                   pi3k_mtor)

hallmark_res <- map_df(res.list,
                   Count_gsea,
                    hallmark)


setdiff(gene.important,unique(mtor_res$gene))
#"CAMK4"  "MKNK1"  "PIK3C3" "PRKCB" These four genes have no effect on the pathway
## To verify the effect of nonexistent genes on the mtor pathway, it could not be verified
hallmark_res.add <- map_df(res.list[18],
                           Count_gsea,
                           hallmark)

  
## View p value
gsea_res$pvalue %>% 
  p.adjust(method = "fdr") %>% 
  p.adjust(method = "fdr") ->p

## Filter interest path
hallmark_res %>% 
  filter(ID=="HALLMARK_PI3K_AKT_MTOR_SIGNALING")->
  mtor_res

hallmark_res %>% 
  filter(ID=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")->
  emt_res

## drawing
mtor_res %>% 
  filter(abs(NES)>1&pvalue<0.05&qvalues<0.25) %>% 
  mutate(group=ifelse(NES>1,"activated","suppressed")) %>% 
  ggplot(aes(reorder(gene,-NES),NES,fill=qvalues))+
  geom_bar(stat = "identity")+
  coord_flip()+
  facet_grid(group~.,scales = "free")+
  scale_fill_gradient(low="#0099F7",high="#f11712")+
  labs(title = "Hallmark PI3K AKT MTOR Signaling",
       x="Gene",
       y="Normalized Enrichment Score")+
  theme_bw(base_size = 18)+
  theme(plot.title = element_text(hjust = 0.5,size=16))
  x11()

# hallmark result drawing ------------------------------------------------------------
  
  gene.important.tab %>% 
    mutate(ind=case_when(HR>1&index=="OS"~"OS risk gene",
                         HR<1&index=="OS"~"OS protective gene",
                         HR>1&index=="DFS"~"DFS risk gene",
                         HR<1&index=="DFS"~"DFS protective gene")) %>% 
    mutate(ind=case_when(gene=="RAC1"~"OS and DFS risk gene",
                         gene%in%c("ADCY2","PIK3R3")~"OS and DFS protective gene",
                         T~as.character(ind)))->
    gene.important.tab.group
  

  
  
    
library(pheatmap)
  
load("dat_output/Results of survival analysis and expression profile analysis of pathway genes.Rdata")

hallmark_res %>%
  pull(ID) %>% 
  table() %>% 
  sort(decreasing = T) %>% 
  as.data.frame() %>% 
  dplyr::rename("pathway"=1)->
  hallmarks.tab



hallmarks.tab %>% 
  filter(Freq>=17) %>%  #The frequency corresponding to the path in bit 5
  pull(pathway) %>% 
  unique()->
  pathway_top5

hallmark_res %>% 
  filter(abs(NES)>1&pvalue<0.05&qvalues<0.25) %>% 
  filter(ID%in%pathway_top5)%>% 
  dplyr::select(gene,ID,NES) %>%
  pivot_wider(names_from = gene,
              values_from = NES) %>% 
  column_to_rownames("ID") %>%  
  replace_na(list(.=0))->
  heat_dat_hall


heat_dat_hall[is.na(heat_dat_hall)] <- 0
heat_dat_hall
names(heat_dat_hall) %>% 
  as.data.frame() %>% 
  set_names("gene") %>% 
  inner_join(gene.important.tab.group %>% 
              dplyr::select(gene,compare,ind)) %>% 
  distinct()->
  anno_col

# Annotated data
anno_col_rev <- data.frame(
  survival_impact=anno_col$ind[seq(1,63,3)],
  paracancer_vs_cancer=anno_col$compare[seq(1,63,3)],
  paracancer_vs_normal=anno_col$compare[seq(2,63,3)],
  cancer_vs_normal=anno_col$compare[seq(3,63,3)]
)

rownames(anno_col_rev) <- names(heat_dat_hall)
rownames(heat_dat_hall) %>% 
  str_remove_all("HALLMARK_") %>% 
  str_replace_all("_"," ") %>% 
  str_to_title() %>% 
  str_replace("Tgf","TGF") %>% 
  str_replace("Il6","IL16") %>% 
  str_replace("Uv","UV") %>% 
  str_replace("Kras","KRAS") %>% 
  str_replace("Dna","DNA") %>% 
  str_replace("Dn","DN") %>% 
  str_replace("Pi3k","PI3K") %>% 
  str_replace("G2m","G2M") %>% 
  str_replace("E2f","E2F") %>% 
  str_replace("Tnfa","TNFA") %>% 
  str_replace("Nfkb","NFKB")->
  rownames(heat_dat_hall)
heat_dat_hall[heat_dat_hall==0] <- NA
 pheatmap(heat_dat_hall,
         annotation_col = anno_col_rev,
         annotation_legend = T,
         cutree_col=2,
         cutree_rows = 3,
         gaps_row = c(10,15,18),
         cellheight = 15,
         cellwidth = 20,
         fontsize=12,
         na_col = "#ffffff",
         show_colnames  = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, 
                                  name ="OrRd")))(100)) ->hp

x11()
library(ggplotify)
ggsave(as.ggplot(hp),filename = "testhm.pdf",width =  13.3 ,height = 7.9)
geneorder <- hp$tree_col$labels[hp$tree_col$order]

# KEGG --------------------------------------------------------------------
kegg <- read.gmt("dat_input/c2.cp.kegg.v7.5.1.symbols.gmt")
kegg_res <- map_df(res.list,
                   Count_gsea,
                   kegg)

kegg_res %>%
  filter(abs(NES)>1&pvalue<0.05&qvalues<0.25) %>% 
  filter(str_detect(ID,"PATHWAY")) %>% 
  pull(ID) %>% 
  table() %>% 
  sort(decreasing = T) %>% 
  as.data.frame() %>% 
  dplyr::rename("pathway"=1)->
  kegg.tab

kegg.tab$Freq %>% 
  unique() %>% 
  .[5]->
  index
kegg.tab$pathway %>% n_distinct() 
kegg.tab %>% 
  filter(Freq>= index) %>%  #The frequency corresponding to the path in bit 5
  pull(pathway) %>% 
  unique()->
  pathway_top5

kegg_res %>% 
  filter(ID%in%pathway_top5) %>% 
  dplyr::select(gene,ID,NES) %>%
  pivot_wider(names_from = gene,
              values_from = NES,
              values_fill = 0) %>% 
  column_to_rownames("ID") %>% 
  dplyr::select(geneorder)->
  heat_dat_kegg


rownames(heat_dat_kegg) %>% 
  str_remove_all("KEGG_") %>% 
  str_replace_all("_"," ") %>% 
  str_to_title() %>% 
  str_replace("Mapk","MAPK") %>% 
  str_replace("Gnrh","GNRH") %>% 
  str_replace('Ppar',"PPAR") %>% 
  str_replace("Jak Stat","JAK STAT")->
  rownames(heat_dat_kegg)

bk <- c(seq(-9,-0.1,by=0.01),seq(0,9,by=0.01))
heat_dat_kegg[heat_dat_kegg==0] <- NA
pheatmap(heat_dat_kegg[geneorder],
         #annotation_col = anno_col_rev[geneorder,],
         #annotation_col=NA,
         annotation_legend = F,
         cluster_cols = F,
         cutree_col=2,
         cutree_rows = 2,
         gaps_row = c(2),
         gaps_col = 11,
         cellheight = 15,
         cellwidth = 20,
         fontsize=12,
         legend_labels = F,
         legend = T,
         show_colnames = F,
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         na_col = "#ffffff")->hp2

ggsave(as.ggplot(hp2),filename = "kegg_nes.pdf",width =  13.3 ,height = 7.9)
x11()




# GO ----------------------------------------------------------------------
go <- read.gmt("dat_input/c5.go.bp.v7.5.1.symbols.gmt")

go_res <- map_df(res.list,
                   Count_gsea,
                 go)
go_res 



go_res %>%
  filter(abs(NES)>1&pvalue<0.05&qvalues<0.25) %>% 
  filter(str_detect(ID,regex("phosphorylation|pi3k",ignore_case = T))) %>% 
  pull(ID) %>% 
  table() %>% 
  sort(decreasing = T) %>% 
  as.data.frame() %>% 
  dplyr::rename("pathway"=1)->
  go.tab

go.tab$Freq %>% 
  unique() %>% 
  .[5]->
  index

go.tab %>% 
  filter(Freq>= index) %>%  #The frequency corresponding to the path in bit 5
  pull(pathway) %>% 
  unique()->
  pathway_top5

go_res %>% 
  filter(ID%in%pathway_top5) %>% 
  dplyr::select(gene,ID,NES) %>%
  pivot_wider(names_from = gene,
              values_from = NES,
              values_fill = 0) %>% 
  column_to_rownames("ID") %>% 
  dplyr::select(geneorder)->
  heat_dat_go


rownames( heat_dat_go) %>% 
  str_remove_all("GOBP_") %>% 
  str_replace_all("_"," ") %>% 
  str_to_title() %>% 
  str_replace("Smad","SMAD") %>% 
  str_replace("Gnrh","GNRH") %>% 
  str_replace('Ppar',"PPAR") %>% 
  str_replace("Jak Stat","JAK STAT")->
  rownames(heat_dat_go)
bk <- c(seq(-9,-0.1,by=0.01),seq(0,9,by=0.01))
heat_dat_go[ heat_dat_go==0] <- NA
pheatmap( heat_dat_go[geneorder],
        # annotation_col = anno_col_rev[geneorder,],
         annotation_legend = F,
         cutree_col=2,
         cutree_rows = 2,
         gaps_row = c(8),
         gaps_col = 11,
         cluster_cols = F,
         cellheight = 15,
         cellwidth = 20,
         fontsize=12,
         color = c(colorRampPalette(colors = c("#195696","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","#bf3338"))(length(bk)/2)),
        na_col = "#FFFFFF")->hp3
x11()

ggsave(as.ggplot(hp3),filename = "go_nes.pdf",width =  13.3 ,height = 7.9)

# immune ------------------------------------------------------------------

immune <- read.gmt("dat_input/c7.all.v7.5.1.symbols.gmt")

immune_res <- map_df(res.list,Count_gsea,
                 immune)



# Write a function to extract heat map data -------------------------------------------------------------
creat.hd <- function(res,topn,geneorder=everything()){
  res %>%
    filter(abs(NES)>1&pvalue<0.05&qvalues<0.25) %>% 
    pull(ID) %>% 
    table() %>% 
    sort(decreasing = T) %>% 
    as.data.frame() %>% 
    dplyr::rename("pathway"=1)->
    res.tab
  
  res.tab$Freq %>% 
    unique() %>% 
    .[topn]->
    index
  
  res.tab %>% 
    filter(Freq>= index) %>%  #The frequency corresponding to the path in bit 5
    pull(pathway) %>% 
    unique()->
    pathway_top5
  
res %>% 
    filter(ID%in%pathway_top5) %>% 
    dplyr::select(gene,ID,NES) %>%
    pivot_wider(names_from = gene,
                values_from = NES,
                values_fill = 0) %>% 
    column_to_rownames("ID") %>% 
    dplyr::select(geneorder)->
    heat_dat
 

heat_dat
  
}



immueheat <- creat.hd(immune_res,5)

# Save the result --------------------------------------------------------------------
save(list=ls(pattern = "res")[-5],file = "dat_output/GSEA enrichment results.Rdata")
load("dat_output/GSEA enrichment results.Rdata")

