library(tidyverse)
library(ranger)
library(randomForest)

load("dat_output/随机森林准备数据.Rdata")
load("dat_output/通路基因生存分析及表达谱分析结果.Rdata")


gene.important.tab %>% 
  filter(index=="OS"&HR>1) %>% 
  filter(gene!="RPS6KA3") %>% 
  pull(gene) %>% 
  unique()->
  os.risk_gene


gene.important.tab %>% 
  filter(index=="DFS"&HR>1) %>% 
  filter(gene!="RPS6KA3") %>% 
  pull(gene) %>% 
  unique()->
  dfs.risk_gene
c( os.risk_gene,dfs.risk_gene) %>% n_distinct()

intersect(os.risk_gene,dfs.risk_gene)

gene.important.tab %>% 
  filter(index=="OS"&HR<1) %>% 
  filter(gene!="RPS6KA3") %>% 
  pull(gene) %>% 
  unique()->
  os.pro_gene


gene.important.tab %>% 
  filter(index=="DFS"&HR<1) %>% 
  filter(gene!="RPS6KA3") %>% 
  pull(gene) %>% 
  unique()->
  dfs.pro_gene
c(os.pro_gene,dfs.pro_gene) %>% n_distinct()

intersect(os.pro_gene,dfs.pro_gene)
# OS case variable ------------------------------------------------------------------

gene.sel <-gene.important.tab %>% 
  filter(index=="OS") %>% 
  filter(gene!="RPS6KA3") %>% 
  pull(gene) %>% 
  unique()

dt.rf <- os.data %>% 
  select(surv_status,gene.sel) %>% 
  rename("OS"=surv_status)
  
  
#Binary classification tree, outcome considering overall survival

ntree <- 500
mtry <- floor(sqrt(length(gene.sel)))
weight <- 0.999999 # The algorithm requires that 1 cannot be taken but is infinitely close to 1. This parameter represents the probability of the variable being selected.
seed <- 001

set.seed(seed) # Setting external seeds ensures repeatable results
surv.rf <- ranger(formula = OS ~ ., 
                  data = dt.rf,
                  num.trees = ntree,
                  mtry = mtry,
                  importance = "impurity",
                  split.select.weights = rep(weight,length(gene.sel)))

# Variables are listed in descending order of importance
var.imp <- sort(ranger::importance(surv.rf),decreasing = T)

var.imp %>% 
  as.data.frame() %>% 
  rownames_to_column("os.gene") %>% 
  rename("vis"=2)->
  os.gene.vis



# dfs ---------------------------------------------------------------------


gene.sel <-gene.important.tab %>% 
  filter(index=="DFS") %>% 
  filter(gene!="RPS6KA3") %>% 
  pull(gene) %>% 
  unique()

dt.rf <- os.data %>% 
  select(surv_status,gene.sel) %>% 
  rename("DFS"=surv_status)


#Binary classification tree, outcome considering overall survival

ntree <- 500
mtry <- floor(sqrt(length(gene.sel)))
weight <- 0.999999 # The algorithm requires that 1 cannot be taken but is infinitely close to 1. This parameter represents the probability of the variable being selected.
seed <- 001

set.seed(seed) # Setting external seeds ensures repeatable results
surv.rf <- ranger(formula = DFS ~ ., 
                  data = dt.rf,
                  num.trees = ntree,
                  mtry = mtry,
                  importance = "impurity",
                  split.select.weights = rep(weight,length(gene.sel)))

# Variables are listed in descending order of importance
var.imp <- sort(ranger::importance(surv.rf),
                decreasing = T)

var.imp %>% 
  as.data.frame() %>% 
  rownames_to_column("dfs.gene") %>% 
  rename("vis"=2)->
  dfs.gene.vis


os.gene.vis %>% 
  mutate(index="os") %>% 
  rename("gene"=os.gene) %>% 
  bind_rows(dfs.gene.vis %>% 
              mutate(index="dfs") %>% 
              rename("gene"=dfs.gene))->
    ranger.res




ranger.res$gene %>% n_distinct()


setwd("C:/Users/Cloudy/Desktop/论文写作/文章写作预分析/table")
gene.important %>% 
  as.data.frame() %>% 
  set_names("gene.important") %>% 
  filter(gene.important!="RPS6KA3") %>% 
  openxlsx::write.xlsx("gene.important.xlsx")

save(os.gene.vis,dfs.gene.vis,ranger.res,file = "dat_output/随机森林算法结果.Rdata")

ranger.res %>% 
ggplot( aes(x= vis, y = gene)) + 
  geom_segment(aes(yend = gene),xend = 0, colour= "grey50",
               size=0.8)+ #The geom_segment function replaces gridlines with "line segments with data points as endpoints.
  geom_point(size = 3, aes(colour = index))+ #Map the lg variable to the color property of the point
  scale_colour_brewer(palette = "Set1")+
  theme_bw(base_size = 15)+
  theme(
    panel.grid.major.y = element_blank(),#Delete the horizontal gridlines
    )+
  facet_grid(index~., scales = "free")





ranger.res %>% 
  filter(index=="os") %>% 
  mutate(group=case_when(gene%in%c(gene_danger.all)~"OS risk gene",
                         T~"OS protective gene")) %>% 
  mutate(group=factor(group,levels = c("OS risk gene","OS protective gene"))) %>% 
  ggplot(aes(x= reorder(gene,vis), y = vis,fill=group)) + 
  geom_bar(stat = "identity")+
  coord_flip()+
  theme_bw(base_size = 15)+
  scale_y_continuous(expand = c(0,0.2))+
  theme_classic(base_size = 18)+
  theme(
    panel.grid.major.y = element_blank(),#Delete the horizontal gridlines
  )+labs(x="Gene",y="Vairable important score")+
  facet_grid(group~.,scales = "free")->a




ranger.res %>% 
  filter(index=="dfs") %>% 
  mutate(group=case_when(gene%in%c(gene_danger.all)~"DFS risk gene",
                         T~"DFS protective gene")) %>% 
  mutate(group=factor(group,levels = c("DFS risk gene","DFS protective gene"))) %>% 
  ggplot(aes(x= reorder(gene,vis), y = vis,fill=group)) + 
  geom_bar(stat = "identity",width = 0.4)+
  scale_fill_brewer(palette = "Set1")+
  coord_flip()+
  theme_bw(base_size = 15)+
  scale_y_continuous(expand = c(0,0.2))+
  theme_classic(base_size = 18)+
  theme(
    panel.grid.major.y = element_blank(),#Delete the horizontal gridlines
  )+labs(x="Gene",y="Vairable important score")->b


library(patchwork)
a+b+
  plot_annotation( tag_levels = "A",
                   title = )+
  plot_layout(guides = 'collect',
              design = "AAB")->c





x11()
