# Correlation of gene expression level
load("dat_output/通路基因生存分析及表达谱分析结果.Rdata")
load("dat_input/NSCLC-exp-clin-data.Rdata")

library(tidyverse)
library(ggcorrplot)
library(ggthemes)
NSCLC_fpkm <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)




gene.important.tab%>% 
  filter(HR>1) %>% 
  pull(gene) %>% 
  unique()->gene.danger


gene.important.tab%>% 
  filter(HR<1) %>% 
  pull(gene) %>% 
  unique()->
  gene.pro

gene.pro%in%gene.danger  #RPS6KA3 results contradict each other and need to be removed.


NSCLC_fpkm %>%
  t() %>% 
  as.data.frame() %>% 
  select(gene.important) %>% 
  select(-"RPS6KA3")->
  gene.imp.exp

cor.res <- cor(gene.imp.exp)

#p value
p.mat <- round(cor_pmat(gene.imp.exp),3)


# As long as the correlation coefficient is greater than or equal to 0.3 above the display
if (F) {
  
  
  # sigmat <- apply(cor.res,2,function(x){res <- ifelse(abs(x)<0.3,T,F) })
  # 
  # p.mat[sigmat] <- 1# Arbitrary assignment，让它大于0.05
}


#Adjust the order so that risk genes come first and protective genes come later
cor.res %>% 
  as.data.frame() %>% 
  select(sort(setdiff(gene.danger,"RPS6KA3")),everything()) %>% 
  .[c(sort(setdiff(gene.danger,"RPS6KA3")),sort(setdiff(gene.pro,"RPS6KA3"))),] %>% 
  as.matrix()->
  cor.res_modi

p.mat %>% 
  as.data.frame() %>% 
  select(sort(setdiff(gene.danger,"RPS6KA3")),everything()) %>% 
  .[c(sort(setdiff(gene.danger,"RPS6KA3")),sort(setdiff(gene.pro,"RPS6KA3"))),] %>% 
  as.matrix()->
  p.mat_modi
  

x11()

ggcorrplot(cor.res_modi,
           method = "circle",
           hc.order = F,
           hc.method = "ward.D",
           outline.color = "white",
           ggtheme =theme_bw(),
           type="lower",
           colors = c("#1e90ff",
                      "white",
                      "#f14040"),
           lab = T,
           lab_size = 3,
           p.mat =  p.mat_modi,
           insig = "blank",
           pch.cex=30,
           tl.cex=rel(1.8))


