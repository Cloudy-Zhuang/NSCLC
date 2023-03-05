
library(tidyverse)
library(data.table)
library(openxlsx)

# miRNA prediction ---------------------------------------------------------------

# Targetscan --------------------------------------------------------------


# The lower the context++score score, the higher the probability that the site is the target; in addition, percentile is the conversion of score, the closer the value is to 100, the higher the probability that the site is the real target。
gene.important <- c("ADCY2", "ARHGAP11A", "CAMK4", "CCNA2", "CCNB2", "CDCA5", "CENPH", "CFL1", "DSCC1", "FASLG", "MAD2L1", "MAPK10", "MKNK1", "PIK3C3", "PIK3R1", "PIK3R3", "PLK4", "PRKCB", "RAC1", "RAF1", "TPX2")

TargetScan_dat <- fread("dat_input/Predicted_Targets_Context_Scores.default_predictions.txt",data.table = F)

TargetScan_dat %>% 
  filter(str_detect(miRNA,"^hsa")) %>% 
  filter(`Gene Symbol`%in%gene.important) %>% 
  filter(`context++ score percentile`>=95) -> # The user-defined parameter is greater than 95
 ts_keep

ts_keep$miRNA %>% n_distinct()

ts_keep %>% 
  group_by(`Gene Symbol`) %>% 
  count() %>% 
  set_names("Target.Gene","n")->ts


# miRTarbase --------------------------------------------------------------
mirtarbase_wr <- read.xlsx("dat_input/miRTarbase/miRTarBase_SE_WR.xlsx")

mirtarbase_wr$Experiments %>% 
  table() %>% 
  as.data.frame() %>% 
  dplyr::rename("type"=1) %>% 
  separate_rows(type,sep="//") %>% 
  filter(type!="") %>% 
  mutate(type=str_remove_all(type,"/")) %>% 
  mutate(type=trimws(type)) %>% 
  distinct(type,.keep_all = T) %>% 
  separate_rows(type,sep = ";") %>% 
  distinct(type,.keep_all = T) %>% 
  arrange(type)->
  type

smirtarbase_wr %>% 
  filter(str_detect(miRNA,"hsa")) %>% 
  filter(Target.Gene%in%gene.important) ->
  mir_tar

mir_tar$miRNA %>% n_distinct()


mirtarbase_wr %>% 
  filter(str_detect(miRNA,"hsa")) %>% 
  filter(Target.Gene%in%gene.important) %>% 
  group_by(Target.Gene) %>% 
  count()->ts2


rbind(ts,ts2) %>% 
pull(Target.Gene) %>% 
  unique()->ind


setdiff(gene.important,ind)



# starbase ----------------------------------------------------------------
setwd("C:/Users/Cloudy/Desktop/论文写作/文章写作预分析/dat_input/starbase")

list.files() %>% 
  str_match(".*_(.*)\\.txt") %>% 
  .[,2]->
  starbase_gene

map_df(list.files(),fread) %>% 
  as.data.frame()->
  starbase_all
starbase_all$clipExpNum %>% summary()

starbase_all %>% 
  filter(clipExpNum>=20) %>% 
  filter(miRanda>=1) %>% 
  distinct(geneName,miRNAname,.keep_all = T)->
  starbase_all_fil

starbase_all_fil$miRNAname %>% n_distinct()

setdiff(gene.important,c(starbase_all_fil$geneName,ind))

c(ts_keep$miRNA,mir_tar$miRNA,starbase_all_fil$miRNAname) %>% n_distinct()



# Build network file ----------------------------------------------------------------
ts_keep %>% 
  select(`Gene Symbol`,miRNA) %>% 
  rename("gene"=`Gene Symbol`) %>% 
  mutate(source="TargetScan") %>% 
  bind_rows(mir_tar %>% 
            select(Target.Gene,miRNA) %>% 
            rename("gene"=Target.Gene) %>% 
            mutate(source="mirtarbase")) %>% 
  bind_rows(starbase_all_fil %>% 
              select(geneName,miRNAname) %>% 
              rename("gene"=geneName,
                     "miRNA"=miRNAname) %>% 
              mutate(source="starbase")) %>% 
  mutate(atr1="miRNA",
         atr2="gene") %>% 
  distinct(gene,miRNA,.keep_all = T)->
  miRNA_net


miRNA_net %>% 
  group_by(miRNA) %>% 
  count() %>% 
  arrange(desc(n))
setwd("C:/Users/Cloudy/Desktop/论文写作/文章写作预分析/table")



miRNA_net %>% 
  write.xlsx("miRNA网络数据.xlsx")
# Transcription factor acquisition ------------------------------------------------------------------


rm(list = ls())
library(httr)
library(jsonlite)

genes = gene.important

url = "https://maayanlab.cloud/chea3/api/enrich/"
encode = "json"
payload = list(query_name = "myQuery", gene_set = genes)

#POST to ChEA3 server
response = POST(url = url, body = payload, encode = encode)
json = content(response, "text")

#results as list of R dataframes
results = fromJSON(json)
ts <- results$`Integrated--meanRank`


results$`Integrated--meanRank` %>% 
  separate_rows(Library,sep = ";") %>% 
  separate_rows(Overlapping_Genes,sep = ",") %>% 
  filter(Overlapping_Genes!="") %>% 
  mutate(Score=as.numeric(Score)) %>% 
  group_split(Overlapping_Genes) %>% 
  map_df(~arrange(.x,Score) %>%
           distinct(TF,.keep_all = T) %>% 
           slice(3)) %>% #slice  just tried out the line numbers, but I can't do anything else
  arrange(Score)  %>% 
  distinct(Overlapping_Genes,TF,.keep_all = T) %>% 
  pull(TF) %>% 
  unique() ->
  TF_top3
 

results$`Integrated--meanRank` %>% 
  separate_rows(Library,sep = ";") %>% 
  separate_rows(Overlapping_Genes,sep = ",") %>% 
  filter(TF%in% TF_top1) ->
  tf_tab



setwd("C:/Users/Cloudy/Desktop/论文写作/文章写作预分析/table")

tf_tab %>% 
  select(TF,Overlapping_Genes) %>% 
  distinct() %>% 
  group_by(TF) %>% 
  count() %>% 
  arrange(desc(n))

tf_tab$TF %>% n_distinct()
  

setwd("D:/桌面重要文件夹/论文写作/文章写作预分析/table")
tf_tab %>% 
  select(TF,Overlapping_Genes) %>% 
  distinct() %>% 
  mutate(atr1="TF",
         atr2="gene") %>% 
  write.xlsx("转录因子网络数据.xlsx")
