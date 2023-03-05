#Consistent clustering
## Clustering algorithm implementation
library(ConsensusClusterPlus)
library(tidyverse)
library(openxlsx)
library(Rtsne)

## Defined pathway up-regulation related genes

gene.up <- sort(c("CFL1", "RAC1", "CCNA2", "PLK4", "CDCA5", "MAD2L1", "TPX2", "CCNB2", "ARHGAP11A", "DSCC1", "CENPH", "FASLG"))

## Pathway downregulates related genes
gene.down <- sort(c("PIK3R3", "ADCY2", "PIK3R1", "MAPK10", "RAF1"))

gene.all <- c(gene.up,gene.down)

## data
load("dat_input/NSCLC-exp-clin-data.Rdata")

##Gene table matrix
NSCLC_fpkm <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)

##Patient ID
datClean <- function(fpkmanno,clin_tab,mygene){
  
  mygene <- mygene[which(mygene%in%rownames(fpkmanno))] 
  
  
  fpkmanno %>% 
    t() %>% 
    as.data.frame() %>% 
    dplyr::mutate(id=substr(rownames(.),1,12)) %>% 
    dplyr::mutate(ind=rowSums(across(where(is.numeric)))) %>% 
    dplyr::arrange(desc(ind)) %>% 
    dplyr::distinct(id,.keep_all = T) %>% 
    dplyr::select(id,any_of(mygene)) %>% 
    dplyr::inner_join(clin_tab %>% 
                        rownames_to_column("id")) %>% 
    column_to_rownames("id") %>% 
    dplyr::select(any_of(sort(mygene)),everything()) %>% 
    dplyr::rename("surv_time"=OS.time,"surv_status"=OS) %>% 
    dplyr::select(surv_time,surv_status,everything())->
    dat_clean 
  dat_clean
}

NSCLC_datclean <- datClean(NSCLC_fpkm,
                           rbind(LUAD_clin_tab,LUSC_clin_tab),
                           gene.all)




dat <- NSCLC_datclean %>% 
  select(3:19) %>% 
  t() %>% 
  as.matrix()->dat




##Normalized to the median, this step is recommended by the official documentation, but some articles do not
dat= sweep(dat,1, apply(dat,1,median,na.rm=T))

results <- ConsensusClusterPlus(dat, 
                                maxK =3, #Maximum number of clusters
                                reps = 1000, # The number of resamples
                                pItem = 0.8,# The resampling ratio of the sample
                                pFeature = 0.8,  # Official default
                                clusterAlg = "km", # The clustering algorithm used can be "hc"(hclust), "pam", "km"(k-means).

                                distance = "euclidean",# The method of distance calculation canbe pearson, spearman, euclidean, binary, maximum, canberra, minkowski
                               # title = title,
                                plot = "png",
                                writeTable=TRUE,
                                seed=123,# Set random seeds for easy repetition
                                innerLinkage="complete",
                                #innerLinkage="average" # The export type of the resulting image can be "png" or "pdf".

)  


# PAC method determines the optimal clustering --------------------------------------------------------------

Kvec = 2:6 #maxK
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}#end for i
# The optimal K
optK = Kvec[which.min(PAC)]

# Clustering into 2 categories is recommended

## Extract data and redraw (code can be adjusted according to the actual situation)


annCol <- data.frame(results = paste0("Cluster",
                                      results[[2]][["consensusClass"]]), #Extract the corresponding category data
                     row.names = colnames(dat))
head(annCol)
library(pheatmap)
library(RColorBrewer)
mycol <- brewer.pal(9, "Set3")
### Can be appended
annColors <- list(Cluster = c("Cluster1" = mycol[1],
                              "Cluster2" = mycol[2]))


annColors <- list(Cluster = c("Cluster1" = mycol[1],
                              "Cluster2" = mycol[2]))

heatdata <- results[[2]][["consensusMatrix"]]
dimnames(heatdata) <- list(colnames(dat),
                           colnames(dat))


pheatmap(mat = heatdata,
         color = colorRampPalette((c("white", "steelblue")))(100),
         border_color = NA,
         annotation_col = annCol,
         annotation_colors = annColors,
         show_colnames = F,
         show_rownames = F)

#To extract the typing result

# The typing results were taken out
Cluster_res <- annCol %>%
  as.data.frame() %>%
  rownames_to_column("patient_ID")



Cluster_res$results %>% table()

##The effect of different categories on prognosis

NSCLC_datclean %>% 
  select(1:2) %>% 
  rownames_to_column("id") %>% 
  inner_join(Cluster_res %>% 
               rename("id"=1)) %>% 
  column_to_rownames("id")->
  NSCLC_cluster


library(survival)
library(survminer)

NSCLC_cluster$group=NSCLC_cluster$results
  
  # Survival analysis 
  # survfit object
  surv.obj = as.formula(paste0("Surv(surv_time, surv_status)~","group"))
  # K-M survival curves
  fit = surv_fit(surv.obj, data = NSCLC_cluster) # 用survival::survfit报错: https://github.com/kassambara/survminer/issues/403
  # log-rank p http://rpkgs.datanovia.com/survminer/reference/surv_pvalue.html
  pval = surv_pvalue(fit)$pval
  # Cox proportional hazards regression model
  cox = coxph(surv.obj, data = NSCLC_cluster)
  cox_summary = summary(cox)
  cox_res = c(cox_summary$coefficients[,c(1:2,5)],
              cox_summary$conf.int[,c(3,4)])


#p values don't make sense

# Give up the plan