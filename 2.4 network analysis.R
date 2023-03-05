## Network analysis

library(tidyverse)
library(openxlsx)
library(data.table)
## Defined pathway up-regulation related genes
gene.up <- sort(c("CFL1", "RAC1", "CCNA2", "PLK4", "CDCA5", "MAD2L1", "TPX2", "CCNB2", "ARHGAP11A", "DSCC1", "CENPH", "FASLG"))


## Pathway downregulates related genes
gene.down <- sort(c("PIK3R3", "ADCY2", "PIK3R1", "MAPK10", "RAF1"))

gene.all <- c(gene.up, gene.down) %>% sort()
## miRNATarbase
mirna_tarbase <- dir(
  path = "dat_input/miRTarbase/",
  pattern = "WR.xlsx$",
  full.names = T,
  recursive = T
)

map_df(mirna_tarbase, read.xlsx) %>% 
  filter(Target.Gene %in% gene.all) %>%
  select(miRNA, Target.Gene) %>% 
  rename("gene"=2)%>%
  mutate(source = "miRNATarbase") ->
miRNATarbase_dat
miRNATarbase_dat$miRNA %>% n_distinct()
miRNATarbase_dat$gene %>% n_distinct()
## starbase
##https://starbase.sysu.edu.cn/tutorial.php
starbase_file <- dir(
  path = "dat_input/starbase/",
  pattern = "txt$",
  full.names = T,
  recursive = T
)

map_df(starbase_file ,fread) %>% 
  filter(clipExpNum>=5&miRanda>=1) %>%   #strict stringency
  select(miRNAname,geneName) %>% 
  filter(geneName%in%gene.all) %>% 
  set_names(c("miRNA","gene")) %>% 
  distinct() %>% 
  mutate(source="starbase")->
  starbase_dat
starbase_dat$miRNA %>% n_distinct()
starbase_dat$gene %>% n_distinct()

## targetscan

starbase_file <- dir(
  path = "dat_input/",
  pattern = "txt$",
  full.names = T,
  recursive = T
)

## target scan
#https://ibook.antpedia.com/x/587474.html#:~:text=%E5%85%B6%E4%B8%AD%EF%BC%8CContext%2B%2B,score%E6%98%AF%E9%A2%84%E6%B5%8B%E7%9A%84%E9%9D%B6%E7%82%B9%E7%9A%84%E7%BB%BC%E5%90%88%E5%88%86%E6%95%B0%EF%BC%8C%E5%88%86%E6%95%B0%E8%B6%8A%E4%BD%8E%E8%AF%81%E6%98%8E%E9%A2%84%E6%B5%8B%E4%BD%8D%E7%82%B9%E4%B8%BA%E7%9C%9F%E5%AE%9E%E9%9D%B6%E7%82%B9%E7%9A%84%E6%A6%82%E7%8E%87%E8%B6%8A%E5%A4%A7%EF%BC%9B%E5%90%8E%E9%9D%A2%E7%9A%84Context%2B%2Bscore%20percentile%E5%88%99%E6%98%AF%E5%89%8D%E9%9D%A2%E5%88%86%E6%95%B0%E7%9A%84%E6%8D%A2%E7%AE%97%EF%BC%8C%E8%B6%8A%E6%8E%A5%E8%BF%91100%EF%BC%8C%E9%A2%84%E6%B5%8B%E4%BD%8D%E7%82%B9%E4%B8%BA%E7%9C%9F%E5%AE%9E%E9%9D%B6%E7%82%B9%E7%9A%84%E6%A6%82%E7%8E%87%E8%B6%8A%E5%A4%A7%E3%80%82
target_scan <- fread("dat_input/Predicted_Targets_Context_Scores.default_predictions.txt",data.table = F)

target_scan %>% 
  filter(`context++ score percentile`>=95) %>% 
  filter(`Gene Symbol`%in%gene.all) %>% 
  filter(str_detect(miRNA,"^hsa")) %>% 
  select(miRNA,`Gene Symbol`) %>% 
  set_names("miRNA","gene") %>% 
  mutate(source="Targetscan")->
  targetscan_dat

targetscan_dat$gene %>% n_distinct()


c(starbase_dat$gene,miRNATarbase_dat$gene,targetscan_dat$gene) %>%
  n_distinct()
gene.all %>% sort()


c(starbase_dat$gene,miRNATarbase_dat$gene,targetscan_dat$gene) %>% 
  unique()->datbase.gene
setdiff(gene.all,datbase.gene)
setdiff(datbase.gene,gene.all)

rbind(miRNATarbase_dat, 
      targetscan_dat,
      starbase_dat) %>% 
  mutate(attr1="miRNA",
         attr2="gene")->
  miRNA_taball
miRNA_taball$miRNA %>% n_distinct()
miRNA_taball %>% 
  group_by(source) %>% 
  summarise(mirna_count=n_distinct(miRNA))

miRNA_taball %>% 
  group_by(miRNA) %>% 
  summarise(gene_count=n_distinct(gene)) %>% 
  arrange(desc(gene_count)) 
  
miRNA_taball$source %>% unique()
miRNA_taball %>% 
 # distinct(miRNA,gene,.keep_all = T) %>%   #Keep it, don't weigh it
  mutate(source_attr1=ifelse(source=="miRNATarbase","Tarbase",NA)) %>% 
  mutate(source_attr2=ifelse(source=="Targetscan","Targetscan",NA)) %>% 
  mutate(source_attr3=ifelse(source=="starbase","starbase",NA))->
  miRNA_taball_unique



# chea3 -------------------------------------------------------------------

rm(list = ls())
library(httr)
library(jsonlite)

genes = gene.all

url = "https://maayanlab.cloud/chea3/api/enrich/"
encode = "json"
payload = list(query_name = "myQuery", gene_set = genes)

#POST to ChEA3 server
response = POST(url = url, body = payload, encode = encode)
json = content(response, "text")

#results as list of R dataframes
results = fromJSON(json)


chea3=results$`Integrated--meanRank`
chea3 %>% 
  separate_rows(Library,sep=";") %>% 
  separate_rows(Overlapping_Genes,sep = ",") %>% 
  filter(Overlapping_Genes!="") %>% 
  mutate(Score=as.numeric(Score)) %>% 
  group_split(Overlapping_Genes) %>% 
  map_df(~arrange(.x,Score) %>% 
           distinct(TF,Overlapping_Genes,.keep_all = T)) %>% 
  group_split(Overlapping_Genes) %>% 
  map_df(~data.frame(Rank=.x$Rank[1:3],TF=.x$TF[1:3],
                     Score=.x$Score[1:3],
                     Overlapping_Genes=.x$Overlapping_Genes[1:3])) %>% 
  mutate(attr1='TF',
         attr2="gene")->
  TF_tab
                    
  
TF_tab %>% 
  group_by(TF) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  DT::datatable()

TF_tab$TF %>% n_distinct()
# write ----------------------------------------------------------------------
write.xlsx(miRNA_taball_unique,"dat_output/miRNA网络数据.xlsx")

write.xlsx(TF_tab,"dat_output/TF网络数据.xlsx")



# cytoscape outputs data -----------------------------------------------------------

TF_node <- read.csv("dat_input/cytoscape网络输出数据/TF node.csv")
library(tidyverse)

dat <- read.xlsx("dat_output/TF网络数据.xlsx")
dat %>% 
  group_by(TF) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  set_names("Transcription factor","Degree") %>%
  remove_rownames() %>% 
  #column_to_rownames("Transcription factor") %>% 
  gridExtra::tableGrob(row=NULL,
                       theme =ggthemes::excel_new_pal())->tab
library(grid)
ttheme_minimal
grid.draw(tab)
x11()



mirna_node <- read.csv("dat_input/cytoscape网络输出数据/Merged Network default.csv")


mirna_node$Degree %>% summary()

mirna_node %>% 
  filter(str_detect(shared.name,"^hsa")) %>% 
  arrange(desc(Degree)) %>% 
  select(shared.name,Degree) %>% 
  slice(1:10) %>% 
  rename("miRNA"=1) %>% 
  gridExtra::tableGrob(row=NULL) %>% 
  grid.draw()
x11()



mirna_node %>% 
  filter(str_detect(shared.name,"^hsa")) %>% 
  arrange(desc(Degree)) %>% 
  pull(shared.name) %>% 
  .[1:3] %>% 
  paste0(collapse = ", ")



hubgene_node <- read.csv("dat_input/cytoscape网络输出数据/string_interactions_short.tsv default.csv")



hubgene_node %>% 
  arrange(desc(Degree)) %>% 
  select(shared.name,Degree) %>% 
  bind_rows(data.frame(shared.name="ADCY2",Degree=0)) %>% 
  rename("Pathway gene"=1) %>% 
  gridExtra::tableGrob(row=NULL) %>% 
  grid.draw()
  x11()

  hubgene_node %>% 
    arrange(desc(Degree)) %>% 
    select(shared.name,Degree) %>% 
    filter(Degree>=7) %>% 
    mutate(ind=paste0(shared.name," (","degree:",Degree,")")) %>% 
    pull(ind) %>% 
    paste(collapse = ", ")
    
