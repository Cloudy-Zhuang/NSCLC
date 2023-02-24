#Function annotation and visualization
##1. KEGG enrichment analysis of PI3K AKT pathway
##2. GSEA analysis of PI3K AKT pathway
##3. GSVA analysis of PI3K AKT pathway

# Load the required package
library(tidyverse)
library(clusterProfiler)
library(GSVA)

R.utils::setOption("clusterProfiler.download.method",'auto')
# 1.KEGG enrichment analysis of PI3K AKT pathway --------------------------------------------------
load("dat_input/TCGA_NSCLC_DEGs.Rdata")


gene <- 
  TCGA_NSCLC_diffgene %>% 
  filter(padj<0.05&abs(log2FoldChange)>=1) %>% 
  pull(gene_id)
#Gene name conversion, return data frame
genename = bitr(gene, 
            fromType="SYMBOL", 
            toType="ENTREZID", 
            OrgDb="org.Hs.eg.db")

EGG <- enrichKEGG(gene = genename$ENTREZID,
                  organism = 'hsa',
                  pvalueCutoff = 0.05)

symboldata <- setReadable(EGG, OrgDb="org.Hs.eg.db", 
                          keyType = "ENTREZID")

symboldata.df <- as.data.frame(symboldata@result)
symboldata.df %>% 
  filter(str_detect(Description,"PI3K"))
symboldata.df %>% 
  filter(str_detect(Description,"MTOR")) # mtor isnt here


## The result was not significant.

## Alternative to the enrich method
load("dat_output/pi3k_inf.Rdata")

pi3k_kegg <- PI3K %>% 
           distinct(genecopy,.keep_all = T) %>% 
           select(pathway.name,
                    gene) %>% 
             rename("term"=1)


kegg_res <- enricher(gene,
                     TERM2GENE = pi3k_kegg,
                     minGSSize = 1,
                     maxGSSize = 20000)
res <- kegg_res@result

## Reads into official files
pi3k_kegg <- read.gmt("dat_input/WP_PI3KAKT_SIGNALING_PATHWAY.v7.5.1.gmt")
## Re-enrichment
kegg_res <- enricher(gene,
                     TERM2GENE = pi3k_kegg,
                     minGSSize = 1,
                     maxGSSize = 20000)
res <- kegg_res@result


## Again the result was not significant and this practice was abandoned.



# ## 2.PI3K Akt pathway GSEA analysis --------------------------------------------------

TCGA_NSCLC_diffgene %>% 
    mutate(rank = rank(log2FoldChange,  
                       ties.method = "random")) %>%
  arrange(desc(rank))->TCGA_NSCLC_diffgene


##Genelist Trilogy
##1. Obtaining gene logFC
geneList <- TCGA_NSCLC_diffgene$log2FoldChange

## 2.Nomenclature
names(geneList) = TCGA_NSCLC_diffgene$gene_id
## 3.Sequencing is important
geneList = sort(geneList, decreasing = TRUE)




head(geneList)

# Need the network, everybody will be crowded, but fast
y <- GSEA(geneList,
          TERM2GENE =pi3k_kegg,
          minGSSize = 1,
          maxGSSize = 1000,
          pvalueCutoff=1)
ydf <- y@result
### Mapping
# The results were not statistically significant. Again, this needs to be abandoned.


# ## 3.PI3K AKT Pathway GSVA Analysis --------------------------------------------------
load("dat_input/NSCLC-exp-clin-data.Rdata")

NSCLC_fpkm <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)
#测序数据选择泊松分布，芯片选择高斯分布
mat <- as.matrix(NSCLC_fpkm)
pi3k_list <- as.list(pi3k_kegg)
kegg1 <- gsva(expr=mat, 
              pi3k_list , 
              kcdf="Poisson",
              method = "gsva",
              parallel.sz=10)

kegg1 %>% 
  t() %>% 
  as.data.frame() %>% 
  rename("pi3k_score"=gene)->
  NSCLC.pi3k.gsva.res
save(NSCLC.pi3k.gsva.res,file="NSCLC_gsva_res.Rdata")


# MOTR pathway------------------------------------------------------------------
## NS

## change to enrich method
mtor <- clusterProfiler::read.gmt("dat_input/HALLMARK_PI3K_AKT_MTOR_SIGNALING.v7.5.1.gmt")




kegg_res <- enricher(gene,
                     TERM2GENE = mtor,
                     minGSSize = 1,
                     maxGSSize = 20000)
res <- kegg_res@result

## NS too



# ## 2.mtor AKT pathway GSEA analysis --------------------------------------------------

TCGA_NSCLC_diffgene %>% 
  mutate(rank = rank(log2FoldChange,  
                     ties.method = "random")) %>%
  arrange(desc(rank))->TCGA_NSCLC_diffgene


## geneList 3steps
## 1.get gene logFC
geneList <- TCGA_NSCLC_diffgene$log2FoldChange

## 2.name
names(geneList) = TCGA_NSCLC_diffgene$gene_id
## 3.sort
geneList = sort(geneList, decreasing = TRUE)




head(geneList)

# need online
y <- GSEA(geneList,
          TERM2GENE=mtor,
          minGSSize = 1,
          maxGSSize = 1000,
          pvalueCutoff=1)
ydf <- y@result
### chart
# NS too. discard.


# ## 3.mtor pathway GSVA analysis --------------------------------------------------
load("dat_input/NSCLC-exp-clin-data.Rdata")

NSCLC_fpkm <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)
# The sequencing data were selected with a Poisson distribution, and the chip was selected with a Gaussian distribution
mat <- as.matrix(NSCLC_fpkm)
mtor_list <- split(mtor$gene, mtor$term)

gsva_res <- gsva(expr=mat, 
                 mtor_list , 
              kcdf="Poisson",
              method = "gsva",
              parallel.sz=10)



