load("dat_output/通路基因生存分析及表达谱分析结果.Rdata")
load("dat_output/Mtor通路研究基因.Rdata")
library(tidyverse)
R.utils::setOption( "clusterProfiler.download.method",'auto' )


library(clusterProfiler)

gene = bitr(gene.important, 
            fromType="SYMBOL", 
            toType="ENTREZID", 
            OrgDb="org.Hs.eg.db")

EGG <- enrichKEGG(gene= gene$ENTREZID,
                  organism     = 'hsa',
                  pvalueCutoff = 1,
                  use_internal_data = T)

EGG <- setReadable(EGG ,
                   OrgDb="org.Hs.eg.db",
                   keyType ="ENTREZID" )

kegg_res <- EGG@result


go <- enrichGO(gene = gene$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all")
go  <- setReadable(go,OrgDb="org.Hs.eg.db",keyType ="ENTREZID" )
go_res <- go@result


go_res %>% 
  group_by(ONTOLOGY) %>% 
  arrange(desc(Count)) %>% 
  dplyr::slice(1:10) %>% 
  dplyr::mutate(generatio=Count/25) %>% 
  mutate(richFactor =Count / as.numeric(sub("/\\d+", "", BgRatio))) %>% 
  arrange(pvalue,desc(richFactor)) %>% 
  ggplot(aes(Description,sort(richFactor),
             color=pvalue,
             size=Count))+
  geom_point()+
  facet_grid(ONTOLOGY~.,scales = "free")+
  coord_flip()+
  scale_color_gradient(low="red",high ="blue" )+
  labs(
    #Set picture title parameters
    x="GO Term",
    y="Rich Factor" ,
    size="Gene Count")+
  theme_bw()+ 
  theme( 
    axis.text.y = element_text(size = rel(1.8)), 
    axis.text.x = element_text(size = rel(1.8)),
    axis.title.x = element_text(size=rel(1.5)),
    axis.title.y = element_text(size = rel(1.5)),
    title = element_text(size = rel(1.5)),
    plot.title = element_text(hjust = 0.5),
    strip.text.y = element_text(size = rel(1.8)))+
  scale_size(range=c(3,8))+
  theme(legend.position = "right")+
  scale_color_viridis_c(begin = 0.6, end = 1)+
  guides(color=guide_colorbar(order=2))->goplot
x11()
print(goplot)

go_res %>% 
  group_by(ONTOLOGY) %>% 
  arrange(desc(Count)) %>% 
  dplyr::slice(1:10) %>% 
  dplyr::select(Description,geneID) %>% 
  gt::gt()


go %>% 
  clusterProfiler.dplyr::group_by(ONTOLOGY) %>% 
  clusterProfiler.dplyr::arrange(desc(Count)) %>% 
  clusterProfiler.dplyr::slice(1:10) %>% 
  barplot()+
  facet_grid(.~ONTOLOGY,scales = "free")+
  scale_fill_continuous(low = "#e06663",
                        high = "#327eba", 
                        name = "pvalue",
                        guide = guide_colorbar(reverse = TRUE, order=1), 
                        trans='log10')



# plot --------------------------------------------------------------------

x = EGG
df = data.frame(x)
## Calculate enrichment fraction
x@result$richFactor =x@result$Count / as.numeric(sub("/\\d+", "", x@result$BgRatio))
y =x@result
showCategory = 20
y %>% 
  filter(pvalue<0.05) %>% 
  nrow()

y %>% 
  filter(pvalue<0.05) %>% 
  dplyr::slice(1:30) %>% 
  #filter(!str_detect(Description,"cancer")) %>% 
  arrange(desc(Count),p.adjust) %>% 
  ggplot(aes(richFactor,forcats::fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=pvalue, size = Count)) +
  ## Adjust the color range. The larger the start, the brighter the overall color
  scale_color_viridis_c(begin = 0.3, end = 1) +
  ## Resize bubbles
  scale_size_continuous(range=c(5, 10)) +
  theme_minimal(base_size = 20) + 
  xlab("Rich factor") +
  ylab("KEGG Pathway") + 
  labs(color="p value")+
  ggtitle("")+
  guides(color=guide_colorbar(order=2))

y %>% 
  filter(pvalue<0.05) %>% 
  filter(!str_detect(Description,"cancer")) %>% 
  arrange(desc(Count),p.adjust) %>% 
  as.data.frame() %>% 
  remove_rownames() %>% 
  dplyr::select(Description,geneID) %>% 
  gt::gt(rownames_to_stub = 1) 



