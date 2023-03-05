## Analysis and visualization of immunoinfiltration results
library(tidyverse)
load("dat_input/nsclc免疫浸润结果.Rdata")
load("dat_input/NSCLC-exp-clin-data.Rdata")
load("dat_input/肺癌1027个癌症样本免疫浸润结果.Rdata")
load("dat_output/TCGA队列通路指数计算结果.Rdata")
load("dat_output/多因素cox分析清洁数据.Rdata")
dat.os.multicox$index_minus %>% 
  na.omit() %>% 
  as.numeric() %>% 
  summary()

dat.os.multicox %>% 
  filter(!is.na(age)) %>% 
  ggplot(aes(age,index_minus,fill=age))+
  geom_boxplot(size=0.75,outlier.colour = NA)+
  ylim(c(-0.3613,1.1))+
  ggpubr::stat_compare_means(size=5,
                        comparisons = list(c("<65","≥65")),
                        label.x = 1.2,
                        label="p.signif",
                      symnum.args =list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                      symbols = c("***", "***", "**", "*", "ns")))+
  scale_fill_manual(values = c("#beb8da","#f6cae5"))+
  theme_bw(base_size = 18)+
  theme(legend.position = "top")+
  labs(y="Pathway index")->age_plot
ggsave(age_plot,filename = "age_pi.pdf",width = 5,height = 5)




## Clinical staging
LUAD_clin <- openxlsx::read.xlsx("dat_input/TCGA-LUAD.GDC_phenotype.xlsx")
LUSC_clin <-data.table::fread("dat_input/TCGA-LUSC.GDC_phenotype.tsv",
                              data.table = F)
LUAD_clin$tumor_stage.diagnoses
LUAD_clin %>% 
  select(submitter_id,tumor_stage.diagnoses) %>% 
  rbind(LUSC_clin %>% 
          select(submitter_id,tumor_stage.diagnoses))->
  stage_info
stage_info$tumor_stage.diagnoses %>% table()
dat.os.multicox %>% 
  rownames_to_column("submitter_id") %>% 
  select(submitter_id,index_minus) %>% 
  inner_join(stage_info) %>% 
  filter(tumor_stage.diagnoses!="not reported") %>% 
  mutate(stage=fct_collapse(tumor_stage.diagnoses,
                            `Stage I`=c("stage i","stage ia",
                                        "stage ib"),
                            `Stage II`=c("stage ii","stage iia",
                                          "stage iib"),
                            `Stage III`=c("stage iii","stage iiia",
                                          "stage iiib"),
                            `Stage IV`=c("stage iv"))) %>% 
  mutate(stage=as.character(stage)) %>% 
  arrange(stage) %>% 
  distinct(submitter_id,.keep_all = T)->stage_dat
stage_dat %>% 
  ggplot(aes(stage,index_minus,fill=stage))+
  geom_boxplot()+
  ggpubr::stat_compare_means(label.y =1.7,
                             label.x = 1,
                             size=4.5)+
  ggpubr::stat_compare_means(comparisons =combn(unique(stage_dat$stage),2,
                                                simplify = F),
                             label = "p.signif")+
  ggsci::scale_fill_jco()+
  theme_bw(base_size = 18)+
  theme(legend.position = "top")+
  labs(x="Stage",y="Pathway index")->a
  
stage_dat %>% 
  mutate(stage=fct_collapse(stage,
                            `Stage I-II`=c("Stage I",
                                           "Stage II"),
                            `Stage III-IV`=c("Stage III",
                                           "Stage IV"))) %>% 
  ggplot(aes(stage,index_minus,fill=stage))+
  geom_boxplot()+
  ggpubr::stat_compare_means(comparisons = list(c("Stage I-II","Stage III-IV")),
                             label = "p.signif")+
  scale_fill_manual(values = c("#fccccc","#f34645"))+
  theme_bw(base_size = 18)+
  theme(legend.position = "top")+
  labs(x="Stage",y="Pathway index")->b

x11()
a
b
library(patchwork)

a+b+patchwork::plot_layout(widths = c(1.25,1))+x11()  


if (F) {
  rt <- dat.os.multicox

  rt$group <- ifelse( rt$index_minus>=0.108,"high","low")
  
  cutoff <- 0.108
  cutnum <- sum(rt$group=="low")
  rt %>% 
    arrange(index_minus) %>% 
    mutate(num=1:nrow(.)) %>% 
    ggplot(aes(x=num,y=index_minus))+
    geom_point(aes(col=group),size=1.5)+
    geom_segment(aes(x = cutnum, y = min(index_minus), 
                     xend = cutnum,
                     yend = cutoff),linetype="dashed")+
    geom_segment(aes(x = 0, y = cutoff, xend = cutnum,yend = cutoff),linetype="dashed")+
    geom_text(aes(x=cutnum/2,y=cutoff+0.15,
                  label=paste0("Cutoff: ",cutoff)),
                  col ="black",size = 4,alpha=0.8)+
    theme_bw(base_size = 18)+
    theme(axis.title.x=element_blank())+
    scale_color_manual(values=rev(c("#2878b5","#ff8884")))+
    scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
    labs(color="group",
         y="Pathway index")->plot.point
  
  ## Plot the distribution of survival events
  
  rt %>% 
    arrange(index_minus) %>% 
    mutate(num=1:nrow(.)) %>%
    mutate(status=case_when(surv_status==0~"Alive",T~"Dead")) %>% 
    ggplot(aes(x=num,y=surv_time))+
    geom_point(aes(col=status),size=1.5)+
    geom_vline(aes(xintercept = cutnum),linetype="dashed")+
    scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
    labs(y="Survival Time (days)",color="OS status")+
    scale_color_manual(values=c("#2878b5","#ff8884"))+
    theme_bw(base_size = 18)+
    theme(axis.title.x=element_blank())->plot.sur
  plot.point/plot.sur->plot_add
  

  x11()
  }






if (T) {
  Compare_means<- function(gene, dat) {
    # Calculate the median value
    median_res <- sort(tapply(dat[,gene],dat$group,median))
    med_compare <-case_when(median_res[1]>median_res[2]~paste0(names(median_res[1])," > ",names(median_res[2])),
                            median_res[1]<median_res[2]~paste0(names(median_res[1])," < ",names(median_res[2])),
                            T~paste0(names(median_res[1]),"=",names(median_res[2]))) 
    
    # Test for homogeneity of variance
    var_check <- var.test(dat[, gene] ~ group, data = dat)
    # If the variance is 0 or cannot be calculated as a missing value, the data will not perform the normal distribution test, so set it to return the missing value directly
    if (var_check$estimate == 0 | is.na(var_check$estimate)) {
      
      return(data.frame(gene = gene, statistic = NA, pvalue = NA, method = NA),
             median.compare=med_compare)
    }
    # Normality test
    shap <- tapply(dat[, gene], dat$group, shapiro.test)
    
    # If any of the two groups of data conform to the normal distribution and the comparison variance of the two groups of data is homogeneous, the T-test is performed; otherwise, wilcox.test is performed
    # Determine whether to perform t test or rank sum test and return statistics and P-values

    if (all(c(shap[[1]]$p.value, shap[[2]]$p.value) > 0.05 && var_check$p.value > 0.05)) {
      test <- t.test(dat[, gene] ~ group, data = dat)
      res <- data.frame(gene = gene, statistic = test$statistic, 
                        pvalue = test$p.value, 
                        method = "t.test",
                        median.compare=med_compare )
    } else {
      test <- wilcox.test(dat[, gene] ~ group, data = dat, exact = F)
      res <- data.frame(gene = gene, statistic = test$statistic, 
                        pvalue = test$p.value, method = "wilcox.test",
                        median.compare=med_compare )
    }
    
    return(res)
  }
  
}



gene.important <- c("ADCY2", "ARHGAP11A", "CAMK4", "CCNA2", "CCNB2", "CDCA5", "CENPH", "CFL1", "DSCC1", "FASLG", "MAD2L1", "MAPK10", "MKNK1", "PIK3C3", "PIK3R1", "PIK3R3", "PLK4", "PRKCB", "RAC1", "RAF1", "TPX2")
risk.gene <- c("ARHGAP11A", "CCNA2", "CCNB2", "CDCA5", "CENPH", "CFL1", "DSCC1", "MAD2L1", "PLK4", "RAC1","TPX2")
pro.gene <-c("ADCY2", "CAMK4", "FASLG", "MAPK10", "MKNK1", "PIK3C3", "PIK3R1", "PIK3R3","PRKCB","RAF1")

pathway.up <- c("CFL1", "RAC1", "CCNA2", "PLK4", "CDCA5", "MAD2L1", "TPX2", "CCNB2", "ARHGAP11A", "DSCC1", "CENPH", "FASLG")
pathway.down <- c("RAF1", "ADCY2", "PIK3R1", "MAPK10", "PIK3R3")

if (F) {
  pathway.reg <- c(pathway.up,pathway.down)
  tab_add <- openxlsx::read.xlsx("1-s2.0-S009286741830237X-mmc1.xlsx",sheet = 2,
                                 startRow = 2)
  tab_add %>% 
    filter(Gene%in%pathway.reg)->
    pathway_regadd
}

setdiff(gene.important,c(risk.gene,pro.gene))


cbind(LUAD_fpkmanno,LUSC_fpkmanno) %>% 
  t() %>% 
  as.data.frame() %>% 
  filter(substr(rownames(.),14,15)!="11") %>% 
  select(sort(pathway.up),sort(pathway.down))->
  gene.exp
  

# ggsea algorithm -----------------------------------------------------------------
genecor <- function(gene,immune_dat){

  immune_dat %>% 
  t() %>% 
  as.data.frame() %>% 
  filter(substr(rownames(.),14,15)!="11") %>% 
   rownames_to_column("id") %>% 
   inner_join(gene.exp %>% 
              rownames_to_column("id") %>% 
                select(id,gene)) %>% 
   select(gene,everything()) %>% 
   column_to_rownames("id")->
   cor_tab

 cor_data <- do.call(rbind,lapply(colnames(cor_tab)[-1],function(x){
   dd <- cor.test(as.numeric(cor_tab[,x]),cor_tab[[gene]],
                  method ="spearman",exact=FALSE)
   data.frame(gene=gene,cell=x,cor=dd$estimate,p.value=dd$p.value)
 }))
 cor_data 

}

cor.res <- map_df(names(gene.exp),
                  genecor,
                  gsva_data) %>% 
  filter(abs(cor)>0.3)

cor.res $pstar <- ifelse(cor.res $p.value < 0.05,
                     ifelse(cor.res $p.value < 0.01,"**","*"),
                     "")

cor.res$gene <- factor(cor.res$gene,
                       levels = names(rev(gene.exp)))
#drawing
#Little details of genetic sequencing
ggplot(cor.res, aes(cell, gene)) + 
  geom_tile(aes(fill = cor), colour = "white",size=0.1)+
  scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  geom_text(aes(label=pstar),col ="black",size = 5)+ 
  theme_minimal(base_size = 18)+# No background
  labs(title = "Correlation of immune infiltration based on ssGSEA")+
  theme(axis.title.x=element_blank(),#No title
        axis.ticks.x=element_blank(),#No X-axis
        axis.title.y=element_blank(),#No y axis
        axis.text.x = element_text(angle = 45, hjust = 1),# Adjust the X-axis text
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=19))+#Adjust the Y-axis text
  #Adjust legend
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))
  
x11()



# cibersort ---------------------------------------------------------------

nsclc_ciber_res %>% 
  select(1:22) %>% 
  t() %>% 
  as.data.frame()->
  ciber_res


ciber_cor.res <- map_df(names(gene.exp),
                        genecor,ciber_res) %>% 
  drop_na() %>% 
  filter(abs(cor)>=0.3)




ciber_cor.res $pstar <- ifelse(ciber_cor.res $p.value < 0.05,
                         ifelse(ciber_cor.res$p.value < 0.01,"**","*"),
                         "")

ciber_cor.res$gene <- factor(ciber_cor.res$gene,
                            levels = names(rev(gene.exp)))
#drawing
ggplot(ciber_cor.res, aes(cell, gene)) + 
  geom_tile(aes(fill = cor), colour = "white",size=0.1)+
  scale_fill_gradient2(low = "#40e0d0",mid = "white",high = "#ff0080")+
  geom_text(aes(label=pstar),col ="black",size = 5)+ 
  theme_minimal(base_size = 18)+# No background
  theme(axis.title.x=element_blank(),#No title
        axis.ticks.x=element_blank(),#No X-axis
        axis.title.y=element_blank(),#No y axis
        axis.text.x = element_text(angle = 45, hjust = 1),# Adjust the X-axis text
        axis.text.y = element_text(size = 12))+#Adjust the Y-axis text
  #Adjust legend
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))+
  labs(title = "Correlation of immune infiltration based on Cibersort")
x11()


# estimate ----------------------------------------------------------------

scores_NSCLC_allcancer %>% 
  mutate(name=str_replace_all(rownames(.),"\\.","-")) %>% 
  remove_rownames() %>% 
  column_to_rownames("name") %>% 
  rownames_to_column("id") %>% 
  inner_join(index_clean %>% 
               select(index_minus) %>% 
               rownames_to_column("id") ) %>% 
  column_to_rownames("id" ) %>% 
  filter(substr(rownames(.),14,15)!="11") %>% 
  mutate(group=ifelse(index_minus>=0.108,"high","low")) -> 
 # select(-index_minus) 
  estimate_index

## Only the patient's estimate index is considered
# estimate_index %>% 
#   rownames_to_column("ID") %>%
#   mutate(ID=substr(ID,1,15)) %>% 
#   distinct(ID,.keep_all = T) %>% 
#   filter(substr(ID,14,15)!="11") %>% 
#   mutate(ID=substr(ID,1,12)) %>% 
#   filter(ID%in%rownames(dat.os.multicox))->
#   estimate_index


cbPalette <- c("#E69F00", "#CC79A7", "#6DBDC4", "#B44541", "#CC79A7", "#F0E442", "#999999","#0072B2","#D55E00")
cols04<-c("#db6968","#4d97cd","#bdc3d2","#6DBDC4",
          "#FA7F6F","#e8c559","#a3d393","#f8984e")
com <- list(c("high","low"))
estimate_index%>%
  filter(substr(rownames(.),14,15)!="11") %>% 
  select(group,StromalScore) %>% 
  ggplot(aes(group,
             StromalScore,
             fill=group,
             color=group))+
  geom_boxplot(alpha=0.2,width=0.3,
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA,
  )+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.75
  )+
  geom_jitter(alpha=0.2,
              position=position_jitterdodge(jitter.width = 0.35, 
                                            jitter.height = 0, 
                                            dodge.width = 0.8))+
   scale_fill_manual(values =
                       cols04[1:2] )+
  scale_color_manual(values = cols04[1:2])+
  ggpubr::stat_compare_means(comparisons = com,
                             label = "p.signif",
                             size=4)+
  # ggpubr::stat_compare_means(size=5,
  #                            label.y = 1.8,
  #                            show.legend = F)+
  theme_bw(base_size = 18)+
  theme(legend.position = "none")+
  #ggsci::scale_fill_jco()+
  #ggsci::scale_color_jco()+
  labs(x="Pathway index",y="Stromal Score")->ss


estimate_index%>%
  select(group,ImmuneScore) %>% 
  ggplot(aes(group,
             ImmuneScore,
             fill=group,
             color=group))+
  geom_boxplot(alpha=0.2,width=0.3,
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA,
  )+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.75
  )+
  geom_jitter(alpha=0.2,
                position=position_jitterdodge(jitter.width = 0.35, 
                                              jitter.height = 0, 
                                              dodge.width = 0.8))+
  scale_fill_manual(values =
                      cols04[3:4] )+
  scale_color_manual(values = cols04[3:4])+
  ggpubr::stat_compare_means(comparisons = com,
                             label = "p.signif",
                             size=4)+
  # ggpubr::stat_compare_means(size=5,
  #                            label.y = 1.8,
  #                            show.legend = F)+
  theme_bw(base_size = 18)+
  theme(legend.position = "none")+
  #ggsci::scale_fill_jco()+
  labs(x="Pathway index",y="Immune Score")->si





estimate_index%>%
  select(group,ESTIMATEScore) %>% 
  ggplot(aes(group,
             ESTIMATEScore,
             fill=group,
             color=group))+
  geom_boxplot(alpha=0.2,width=0.3,
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA,
  )+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.75
  )+
  geom_jitter(alpha=0.2,
              position=position_jitterdodge(jitter.width = 0.35, 
                                            jitter.height = 0, 
                                            dodge.width = 0.8))+
  scale_fill_manual(values =
                      cols04[5:6] )+
  scale_color_manual(values = cols04[5:6])+
  ggpubr::stat_compare_means(comparisons = com,
                             label = "p.signif",
                             size=4)+
  # ggpubr::stat_compare_means(size=5,
  #                            label.y = 1.8,
  #                            show.legend = F)+
  theme_bw(base_size = 18)+
  theme(legend.position = "none")+
  # ggsci::scale_fill_igv()+
  # ggsci::scale_color_igv()+
  labs(x="Pathway index",y="ESTIMATE Score")->se

estimate_index%>%
  select(group,tumor_purity_score) %>% 
  ggplot(aes(group,
             tumor_purity_score,
             fill=group,
             color=group))+
  geom_boxplot(alpha=0.2,width=0.3,
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA,
  )+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.75
  )+
  geom_jitter(alpha=0.2,
              position=position_jitterdodge(jitter.width = 0.35, 
                                            jitter.height = 0, 
                                            dodge.width = 0.8))+
  scale_fill_manual(values =
                      cols04[7:8] )+
  scale_color_manual(values = cols04[7:8])+
  ggpubr::stat_compare_means(comparisons = com,
                             label = "p.signif",
                             size=4)+
  # ggpubr::stat_compare_means(size=5,
  #                            label.y = 1.8,
  #                            show.legend = F)+
  theme_bw(base_size = 18)+
  theme(legend.position = "none")+
  # ggsci::scale_fill_tron()+
  # ggsci::scale_color_tron()+
  labs(x="Pathway index",y="Tumor Purity Score")->sp



library(patchwork)
ss+si+se+sp+
  plot_layout(ncol=4)
x11()
# Relationship between pathway index and immune infiltration ------------------------------------------------------------


indexcor <- function(index,immune_dat){
  
  immune_dat %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("id") %>% 
    inner_join(index_clean %>% 
                 rownames_to_column("id") %>% 
                 select(id,index)) %>% 
    select(index,everything()) %>% 
    column_to_rownames("id")->
    cor_tab
  
  cor_data <- do.call(rbind,lapply(colnames(cor_tab)[-1],function(x){
    dd <- cor.test(as.numeric(cor_tab[,x]),cor_tab[[index]],
                   method ="spearman",exact=FALSE)
    data.frame(index=index,cell=x,cor=dd$estimate,p.value=dd$p.value)
  }))
  cor_data 
  
}




index.ssgsea.res <- indexcor(immune_dat = gsva_data %>% 
                               select(-which(substr(names(.),14,15)=="11")),
                             index = "index_minus")
index.ciber.res <- indexcor(immune_dat = ciber_res %>% 
                           select(-which(substr(names(.),14,15)=="11")),index = "index_minus")


# ## Adjustment of gsva_data to keep only the tumor samples
gsva_data %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  filter(substr(ID,14,15)!="11") %>% 
  # mutate(ID=substr(ID,1,15)) %>%
  # distinct(ID,.keep_all = T) %>%
  # mutate(ID=substr(ID,1,12)) %>%
  # filter(ID%in%rownames(dat.os.multicox)) %>%
  distinct(ID,.keep_all = T) %>%
  column_to_rownames("ID") %>%
  t() %>%
  as.data.frame()->gsva_data_patient



index.ssgsea.res %>% 
  filter(p.value<0.05) %>% 
  #filter(abs(cor)>=0.3) %>% 
  ggplot(aes(cor,forcats::fct_reorder(cell,cor)))+
  geom_segment(aes(xend=0,yend=cell))+
  geom_point(aes(col=p.value,size=abs(cor)))+
  scale_colour_gradientn(colours=c("#108dc7","#ef8e38"))+
  #scale_color_viridis_c(begin = 0.3, end = 1)+
  scale_size_continuous(range =c(3,6))+
  guides(size = guide_legend(order = 1))+
  theme_bw(base_size = 20)+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(size="Cor",
       y="Cell type",
       x="Cor",
       title = "ssGSEA")->
  a



index.ciber.res %>% 
  filter(p.value<0.05) %>% 
  #filter(abs(cor)>=0.3) %>%
  ggplot(aes(cor,forcats::fct_reorder(cell,cor)))+
  geom_segment(aes(xend=0,yend=cell))+
  geom_point(aes(col=p.value,size=abs(cor)))+
  #scale_colour_gradientn(colours=c("#108dc7","#ef8e38"))+
  scale_colour_gradientn(colours=c("#7fc97f","#984ea3"))+
  #scale_color_viridis_c(begin = 0.3, end = 1)+
  scale_size_continuous(range =c(3,6))+
  guides(size = guide_legend(order = 1))+
  theme_bw(base_size = 20)+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(size="Cor",
       y="Cell type",
       x="Cor",
       title = "Cibersort")->
  b

x11()
a

library(patchwork)
x11()
a+b


# Immunoinfiltration results of different pathway index groups ---------------------------------------------------------
library(pheatmap)
index_clean %>% 
  mutate(group=ifelse(index_minus>=0.108,
         "High-Pathway index","Low-Pathway index")) %>% 
  select(group) %>% 
  arrange(group)->anno_col

filter(anno_col,group=="High-Pathway index") %>% 
  rownames()->high_group


pheatmap(select(gsva_data,high_group,everything()),
         annotation_col = anno_col,
         show_colnames = F,
         cluster_cols = F)


pheatmap(select(ciber_res,high_group,everything()),
         annotation_col = anno_col,
         show_colnames = F,
         cluster_cols = F)

ggplot(aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(position = "stack",stat = "identity") + 
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(23))



pheatmap(t(gsva_data),
         show_rownames = F)
# It didn't work out well


# Diagram of different pathway index boxes -------------------------------------------------------

gsva_data %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("id") %>% 
  inner_join(index_clean %>%   #
               select(index_minus) %>% 
               rownames_to_column("id") 
             ) %>% 
  column_to_rownames("id") %>% 
  # rownames_to_column("ID") %>% 
  # mutate(ID=substr(ID,1,15)) %>%
  # distinct(ID,.keep_all = T) %>%
  # mutate(ID=substr(ID,1,12)) %>%
  # filter(ID%in%rownames(dat.os.multicox)) %>%
  # distinct(ID,.keep_all = T) %>%
  # column_to_rownames("ID") %>% 
  mutate(group=ifelse(index_minus>=0.108,"high","low")) %>% 
  select(-index_minus) %>% 
  select(group,everything())->
  gsva_index_patient


cell_com <- function(cell,dat){
    
   dat %>% 
   select(cell,group) %>% 
   Compare_means(cell,.)->res
   
     res
    
}
which(names(ts)=="Natural killer cell")
names(ts)[21] <- "nkcell"
ggpubr::compare_means(nkcell~group,ts)
                     
tapply(ts$nkcell, ts$group, median)


gsva_index_patient %>% 
  filter(substr(rownames(.),14,15)!="11") %>% 
  map_df(names(.)[-1],
         cell_com,
         dat=.) %>% 
  filter(pvalue<0.05&median.compare!="high=low")->
  cell.com.res


cell.com.res$median.compare <- factor(cell.com.res$median.compare,
                                       levels = c("low < high","high < low"))
cell.com.res <- arrange(  cell.com.res,median.compare)



gsva_index_patient %>% 
  filter(substr(rownames(.),14,15)!="11") %>% 
  pivot_longer(names_to = "celltype",
               values_to = "NES",
               -group) %>% 
  filter(celltype%in%cell.com.res$gene) %>% 
  mutate(celltype=factor(celltype,levels = cell.com.res$gene)) %>% 
  ggplot(aes(x = celltype, y = NES))+
  geom_boxplot(aes(fill = group),position = position_dodge(1),width=.3,outlier.shape = NA)+
  geom_violin(aes(colour = group),position = position_dodge(1),scale = "width",fill=NA)+
  theme_bw(base_size = 18)+
  labs(fill="Pathway index",
       color="Pathway index",
       x="Cell type",
       title = "ssGSEA")+
  scale_fill_manual(values = rev(c("#8ecfc9","#ffbe7a")))+
  scale_color_manual(values = rev(c("#8ecfc9","#ffbe7a")))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                   vjust = 1, colour = "black"),
        legend.position = "top")+
  ggpubr::stat_compare_means(aes(group=group), label = "p.signif")->
  c





## extractive
ciber_res %>% 
  t() %>% 
  as.data.frame() %>% 
  #filter(substr(rownames(.),14,15)!="11") %>% 
  rownames_to_column("id") %>% 
  inner_join(index_clean %>%   #
               select(index_minus) %>% 
               rownames_to_column("id") 
  ) %>% 
  column_to_rownames("id") %>% 
  rownames_to_column("ID") %>%
  mutate(ID=substr(ID,1,15)) %>%
  distinct(ID,.keep_all = T) %>%
  mutate(ID=substr(ID,1,12)) %>%
  filter(ID%in%rownames(dat.os.multicox)) %>%
  distinct(ID,.keep_all = T) %>%
  column_to_rownames("ID") %>%
  mutate(group=ifelse(index_minus>=0.108,"high","low")) %>% 
  select(-index_minus) %>% 
  select(group,everything())->
  ciber_index_patient


ciber_index_patient %>% 
  filter(substr(rownames(.),14,15)!="11") %>% 
  map_df(names(.)[-1],
         cell_com,
         dat=.) %>% 
  filter(pvalue<0.05&median.compare!="high=low") %>% 
  arrange(desc(median.compare))->
  ciber.com.res

ciber.com.res$median.compare <- factor(ciber.com.res$median.compare,
                                      levels = c("low < high","high < low"))
ciber.com.res <- arrange(ciber.com.res,median.compare,gene)

ciber_index_patient  %>% 
filter(substr(rownames(.),14,15)!="11") %>% 
pivot_longer(names_to = "celltype",
               values_to = "Proportion",
             -group) %>% 
filter(celltype%in%ciber.com.res$gene) %>% 
mutate(celltype=factor(celltype,levels = ciber.com.res$gene)) %>% 
ggplot(aes(x = celltype, y = Proportion))+
  geom_boxplot(aes(fill = group),position = position_dodge(1),width=.3,outlier.shape = NA)+
  geom_violin(aes(colour = group),position = position_dodge(1),scale = "width",fill=NA)+
  theme_bw(base_size = 18)+
  labs(fill="Pathway index",
       color="Pathway index",
       x="Cell type",
       title = "Cibersort")+
  # scale_fill_manual(values = c("#ffbe7a","#8ecfc9"))+
  # scale_color_manual(values = c("#ffbe7a","#8ecfc9"))+
   scale_fill_manual(values = c("#5fc861","#f9ea28"))+
   scale_color_manual(values = c("#5fc861","#f9ea28"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                   vjust = 1, colour = "black"),
        legend.position = "top")+
  ggpubr::stat_compare_means(aes(group=group), 
                     label = "p.signif")->d


a
b
c
d
x11()

ciber.com.res %>% 
  filter(median.compare=="high < low") %>% 
  pull(gene)->cel1

cell.com.res %>% 
  filter(median.compare=="high < low") %>% 
  pull(gene)->cel2

c(cel1,cel2) %>% sort()

### Bar chart of cibersort infiltration, too many to use
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))

ciber_index_patient %>%
 select(group,  ciber.com.res$gene) %>% 
 filter(group=="high") %>% 
 select(-group) %>% 
 rownames_to_column("sample") %>% 
 pivot_longer(cols =-sample,
              names_to = "celltype",
              values_to = "proportion") %>% 
  ggplot(aes(sample,proportion,fill = celltype)) + 
  geom_bar(position = "stack",stat = "identity") + 
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(23))


# The expression of immune checkpoint genes -------------------------------------------------------------
gene_id <- c("PD-1"="PDCD1",
             "PD-L1"="CD274",
             "HAVCR2"="TIM3",
             "OX40"="TNFRSF4")
gene <- c(
          "PDCD1",
          "CD274",
          "CTLA4",
          "BTLA",
          "TIGIT",
          "HAVCR2",#"TIM3"
          "TNFRSF4" #ox40
) 

cbind(LUAD_fpkmanno,LUSC_fpkmanno) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("id") %>% 
  inner_join(index_clean %>% 
               select(index_minus) %>% 
              rownames_to_column("id") ) %>% 
  column_to_rownames("id" ) %>% 
  filter(substr(rownames(.),14,15)!="11") %>% 
  select(any_of(gene),"index_minus")->ts #Extraction of interest gene



# Tide Database data preparation ------------------------------------------------------------


cbind(LUAD_fpkmanno,LUSC_fpkmanno) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("id") %>% 
  inner_join(index_clean %>% 
               select(index_minus) %>% 
               rownames_to_column("id") ) %>% 
  column_to_rownames("id" ) %>% 
  mutate(group=ifelse(index_minus>=0.108,"high","low")) %>% 
  select(-index_minus) %>% 
  rownames_to_column("id") %>% 
  group_split(group)->
  dat_tide #The rows are the samples and the columns are the genes

# The output file must have rows for genes and columns for samples

Tidepre <- function(dat,filename){
TIDE <- dat
TIDE <- column_to_rownames(TIDE,"id")
TIDE <- as.data.frame(t(TIDE))
TIDE <-  future_apply(TIDE,2,as.numeric)
# # In order to get a better result, two directions are used for median centered

# TIDE <- sweep(TIDE,2, future_apply(TIDE,2,median,na.rm=T))
# TIDE <- sweep(TIDE,1, future_apply(TIDE,1,median,na.rm=T))

write.table(TIDE,filename,sep = "\t",row.names = T,col.names = NA,quote = F)
}

library(future)
library(future.apply)
plan(multisession)
Tidepre( dat_tide[[1]],"index_high_notrev")
Tidepre( dat_tide[[2]],"index_low_notrev")
dat=dat_tide[[1]]


myts <- ts[[1]][1:20]





# immune AI ---------------------------------------------------------------
if(F){
  cbind(LUAD_fpkmanno,LUSC_fpkmanno) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("id") %>% 
    inner_join(index_clean %>% 
                 select(index_minus) %>% 
                 rownames_to_column("id") ) %>% 
    column_to_rownames("id" ) %>% 
    mutate(group=ifelse(index_minus>=0.108,"high","low")) %>% 
    select(-index_minus) %>% 
    rownames_to_column("id") -> 
    dat_ai

  dat_ai %>% 
    mutate(ID=substr(id,1,15)) %>%
    distinct(id,.keep_all = T) %>%
    mutate(id=substr(ID,1,12)) %>%
    filter(id%in%rownames(dat.os.multicox)) %>%
    distinct(id,.keep_all = T) %>%
    column_to_rownames("id") %>% 
    select(group,everything()) %>% 
    t() %>% 
    as.data.frame()->
    dat_ai_patient 

library(ImmuCellAI)

tag <- dat_ai_patient[1,]

# dat_ai_patient %>% 
#   .[-1,] %>% 
#   mutate_all(as.numeric) %>% 
#   mutate_all(~(.x^2-1)) %>% 
#   rbind(tag)->
#   dat_ai_patient_fpmk

dat_ai_patient %>% 
  t() %>% 
  as.data.frame() %>% 
  select(group,everything()) %>%
  select(-ID) %>% 
  t() %>% 
  as.data.frame()->
  dat_ai_patient_rev


dat_ai_patient_rev[1,] <- ifelse(dat_ai_patient_rev[1,]=="high","Group1","Group2")

dat_ai_patient_rev %>% 
rownames_to_column("Symbol") %>% 
data.table::fwrite( "dat_pre.txt")
dat_test <- data.table::fread("dat_pre.txt")



test <- dat_ai_patient_fpmk_rev[1:20]
ts1 <- ImmuCellAI_new(
     tss,
    data_type = "rnaseq",
    group_tag=1,
    response_tag=1
  )

tttt <- data.table::fread("dat_pre.txt",data.table = F)
tss <- data.table::fread("ImmuCellAI_example.txt",data.table = F) %>% 
  column_to_rownames("Symbol")
# Group name replacement
dat_ai_patient[1,] <- ifelse(dat_ai_patient[1,]=="high","Group1","Group2")
dat_ai_patient$`TCGA-05-4250`

ts2 <- ImmuCellAI_new(
  dat_ai_patient[1:10],
  data_type = "ran-seq",
  group_tag=1,
  response_tag=1
)


  }

# ##The significance of correlation analysis is examined -----------------------------------------------------------


library(Hmisc)

res <- rcorr(as.matrix(ts),type="spearman")

CorMatrix <- function(cor,p) {
  ut <- upper.tri(cor)
  data.frame(row = rownames(cor)[row(cor)[ut]] ,
             column = rownames(cor)[col(cor)[ut]], 
             cor =(cor)[ut], 
             p = p[ut] ) } 

cor_dat <- CorMatrix(res$r, res$P) %>% 
  arrange(row)
cor_dat %>%  
  # filter(p<0.05) %>% 
  filter(row=="index_minus"|column=="index_minus")->cor_sig
cor_sig %>% 
  gt::gt()





cor_dat %>%  
  filter(p<0.05) %>% 
  filter(row%in%unique(c(cor_sig$row,cor_sig$column))&column%in%unique(c(cor_sig$row,cor_sig$column)))



unique(c(cor_sig$row,cor_sig$column)) %>% 
  sort() %>% 
  paste0(collapse = ",")

##The genes with statistically different correlation with the target genes were screened
ts %>% 
  select(unique(c(cor_sig$row,cor_sig$column)))->ts_sig

##All reserved presentation
gene_cor <- cor(ts, method = 'spearman')
#Correlation analysis of gene expression values, taking Pearson correlation coefficient as an example
gene_cor <- cor(ts_sig, method = 'spearman')



#Remove the autocorrelation of the gene, which is the value of the diagonal
diag(gene_cor) <- 0
gene_cor  #The final Pearson correlation matrix of inter-gene expression values

#The obtained correlation matrix is converted into pairwise corresponding data frame structure
gene_cor <- reshape2::melt(gene_cor)
gene_cor <- subset(gene_cor, value != 0)  #Remove the correlation of 0 values
head(gene_cor)  #The first two columns are the names of the two genes, and the third is the correlation of the two genes
#chorography
library(circlize)
x11()
par(cex = 2.2, mar = c(0, 0, 0, 0))

chordDiagram(gene_cor, 
             annotationTrack = c('grid', 'name', 'axis'), #Draw the peripheral arc area, showing the name and scale axis
             grid.col = c(CACNA1C = 'red'),
             #grid.col = c(CACNA1C = 'green3', PLVAP = 'red', CDKN3 = 'orange', CDC25C = 'purple', UBE2T = 'skyblue', SKA1 = 'blue'), #Defining gene color
             col = colorRamp2(c(-1, 0, 1), c('green', 'white', 'red'), transparency = 0.5), #Show the color range of the lines according to the size of the correlation
             annotationTrackHeight = c(0.05, 0.05), #Name the distance from the arc, and the width of the arc
)

i <- seq(0,0.995,0.005)
rect(-1+i/2, #xleft
     -1, #ybottom
     -0.9975+i/2, #xright
     -0.96, #ytop
     col = paste(as.character(color[,1]), "FF", sep = ""),
     border = paste(as.character(color[,1]), "FF", sep = ""))
text(-0.97, -1.03, "-1")
text(-0.51, -1.03, "1")
if (F) {
  # Little Y drawing version ------------------------------------------------------------------
  
  Links <- gene_cor %>% 
    arrange(Var1)
  
  
  
  Links$Var1 <- as.character(Links$Var1)
  Links$Var2<- as.character(Links$Var2)
  names(Links) <- c("Gene_1", "Gene_2", "Correlation")
  Corr <- data.frame(rbind(cbind(Links[,1], Links[,3]), cbind(Links[,2], Links[,3])))
  colnames(Corr) <- c("Gene","Correlation")
  
  #Record the original sequence of genes in the Index column
  Corr$Index <- seq(1,nrow(Corr),1)
  #Sort by gene name
  Corr <- Corr[order(Corr[,1]),]
  str(Corr)#View data type
  
  Corr[,2] <- as.numeric(as.character(Corr[,2]))
  
  #Write gene name
  GeneID <- data.frame(GeneID=unique(Corr[,1]))
  
  #There are n plus 1 genes, which are related to n genes

  n <- length(unique(Corr[,1]))-1
  
  #I write the start site of the gene, always starting at 0
  GeneID$Gene_Start <-rep(0,n+1)
  #The sum of correlation coefficients of each gene was calculated in turn as the gene termination site

  for (i in 1:(n+1)){
    GeneID$Gene_End[i] <- sum(abs(Corr[,2])[(i*n - (n - 1)):(i*n)])
  }
  
  #Define enough colors
  mycol <- c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
  #Give each gene a color
  GeneID$Color<-mycol[1:(n+1)]
  head(GeneID)
  
  
  #The width of the line is the absolute value of the correlation coefficient
  Corr[,2] <- abs(Corr[,2]) #Absolute value
  
  for (i in 1 : (n + 1)){
    sub <- as.data.frame(Corr[(i*n - (n - 1)):(i*n),2])
    for (j in 2 : n){
      sub[1,2] <- 0
      sub[j,2] <- sub[j-1,2] + sub[j-1,1]
      sub[1,3] <- sub[1,1]
      sub[j,3] <- sub[j-1,3] + sub[j,1]
    }
    Corr[(i*n - (n - 1)):(i*n),4] <- sub[,2]
    Corr[(i*n - (n - 1)):(i*n),5] <- sub[,3]
  }
  #The genes are returned to their original sequence according to the Index column
  Corr <- Corr[order(Corr$Index),] 
  head(Corr)
  Corr <- drop_na(Corr)
  
  #V4 is the starting position and V5 is the ending position
  #把它写入Links里，start_1和end_1对应Gene_1，start_2和end_2对应Gene_2
  Links$start_1 <- Corr$V4[1:(nrow(Corr)/2)]
  Links$end_1 <- Corr$V5[1:(nrow(Corr)/2)]
  Links$start_2 <- Corr$V4[(nrow(Corr)/2 + 1):nrow(Corr)]
  Links$end_2 <- Corr$V5[(nrow(Corr)/2 + 1):nrow(Corr)]
  head(Links)
  
  #Define the color of the correlation coefficient
  #The maximum correlation coefficient is 1 and the minimum is -1. 201 colors are set here
  #Minus 1 to 0 is the first 100, 0 to 1 is the last 100
  color <- data.frame(colorRampPalette(c("#67BE54", "#FFFFFF", "#F82C2B"))(201))
  #According to the value of the correlation coefficient, the corresponding color is given

  for (i in 1:nrow(Links)){
    Links[i,8] <- substring(color[Links[i,3] * 100 + 101, 1], 1, 7)
  }
  names(Links)[8] <- "color"
  head(Links)
  
  
  if (T) {
    library(circlize)
    
    
    #Plot setting
    pdf("correlation.pdf",width = 6,height = 6)
    
    par(mar=rep(0,4))
    circos.clear()
    circos.par(start.degree = 90, #Where do I start, in counterclockwise order
               #gap.degree specifies the left and right spacing, and track.margin specifies the upper and lower spacing
               gap.degree = 5, #The size of the spacing between gene bars
               track.margin = c(0,0.23), #The larger the value, the smaller the distance between the gene and the line
               cell.padding = c(0,0,0,0)
    )
    
    circos.initialize(factors = GeneID$GeneID,
                      xlim = cbind(GeneID$Gene_Start, GeneID$Gene_End))
    
    #Predraw gene
    circos.trackPlotRegion(ylim = c(0, 1), factors = GeneID$GeneID, 
                           track.height = 0.05, #Genetic lines of fat and thin
                           #panel.fun for each sector
                           panel.fun = function(x, y) {
                             #select details of current sector
                             name = get.cell.meta.data("sector.index") #Gene ID
                             i = get.cell.meta.data("sector.numeric.index") #Gene number
                             xlim = get.cell.meta.data("xlim")
                             ylim = get.cell.meta.data("ylim")
                             
                             #Gene name
                             circos.text(x = mean(xlim), y = 1,
                                         labels = name,
                                         cex = 1.2, #Gene ID text size
                                         niceFacing = TRUE, #Keep the gene name head up
                                         facing = "bending",#Gene names go in an arc direction, which is reverse.clockwise
                                         adj = c(0.5, -2.5), #The location of the gene name controls left and right and up and down, respectively以上翻译结果来自有道神经网络翻译（YNMT）· 通用场景

                                         font = 2  #bold
                             )
                             
                             #plot main sector
                             circos.rect(xleft = xlim[1], 
                                         ybottom = ylim[1],
                                         xright = xlim[2], 
                                         ytop = ylim[2],
                                         col = GeneID$Color[i],
                                         border = GeneID$Color[i])
                             
                             #plot axis
                             circos.axis(labels.cex = 0.7, 
                                         direction = "outside"#,
                                         #The default is nice and you can tweak it with the following parameters
                                         #Master scale setting
                                         #major.at = seq(from = 0,
                                         #               to = floor(GeneID$Gene_End)[i], 
                                         #               by = 400), #It adjusts itself according to the length of the gene

                                         #Subscale quantity
                                         #minor.ticks = 4, 
                                         #labels.niceFacing = TRUE,
                                         #labels.facing = "outside" #or clockwise
                             )})
    
    
    
    #Draw line
    for(i in 1:nrow(Links)){
      circos.link(sector.index1 = Links$Gene_1[i], 
                  point1 = c(Links[i, 4], Links[i, 5]),
                  sector.index2 = Links$Gene_2[i], 
                  point2 = c(Links[i, 6], Links[i, 7]),
                  col = paste(Links$color[i], "C9", sep = ""), 
                  border = FALSE, 
                  rou = 0.7   #links y values for starting and ending points (percentage of circle radius)
      )}
    
    #Drawing example
    i <- seq(0,0.995,0.005)
    rect(-1+i/2, #xleft
         -1, #ybottom
         -0.9975+i/2, #xright
         -0.96, #ytop
         col = paste(as.character(color[,1]), "FF", sep = ""),
         border = paste(as.character(color[,1]), "FF", sep = ""))
    text(-0.97, -1.03, "-1")
    text(-0.51, -1.03, "1")
    dev.off()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  }
  
}
# perforin ---------------------------------------------------------------------
#CYT is defined as the mean value of GZMA and PRF1 expression in TPM
library(patchwork)
# Other immune-blocking related genes
#25 genes
gene_extend <- c("BTN2A2", "IDO1", "TDO2", "VTCN1", "ADORA2A", "C10orf54", "CD276", "CD274","PDCD1", "CTLA4", "CD160", "BTLA", "HAVCR2", "LGALS9", "CD47", "SIRPA", "TIGIT", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR3DL2", "LAG3", "KIR2DL5A", "KIR2DL5B", "KIR3DL1", "KIR3DL3")

cbind(LUAD_fpkmanno,LUSC_fpkmanno) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("id") %>% 
  inner_join(index_clean %>%   #
               select(index_minus) %>% 
               rownames_to_column("id") 
  ) %>% 
  column_to_rownames("id") %>% 
  # rownames_to_column("ID") %>% 
  # mutate(ID=substr(ID,1,15)) %>%
  # mutate(ID=substr(ID,1,12)) %>%
  # filter(ID%in%rownames(dat.os.multicox)) %>%  
  # distinct(ID,.keep_all = T) %>% 
  # column_to_rownames("ID") %>% 
  mutate(group=ifelse(index_minus>=0.108,"high","low")) %>% 
  select(group,everything()) %>% 
  select("GZMA","PRF1",any_of(gene_extend),"index_minus")->
  dat_int_gene







dat_int_gene %>% 
  # rownames_to_column("ID") %>% 
  # mutate(ID=substr(ID,1,15)) %>%
  # mutate(ID=substr(ID,1,12)) %>%
  # filter(ID%in%rownames(dat.os.multicox)) %>%  
  # distinct(ID,.keep_all = T) %>% 
  # column_to_rownames("ID") %>% 
  filter(substr(rownames(.),14,15)!="11") %>% 
  mutate(CYT=(GZMA+PRF1)/2) %>% 
  select(-1,-2) %>% 
  mutate(group=ifelse(index_minus>=0.108,"high","low")) %>% 
  select(-index_minus) %>% 
  select(group,everything())->
  ts.exp

ts.exp %>% 
  ggplot(aes(group,CYT))+
  geom_boxplot()+
  ggpubr::stat_compare_means()

  # Immune block gene expression matrix

setdiff(gene_extend,names(ts.exp))


sort(gene_extend) %>% 
  paste0(collapse = ",")
# Statistical result
com_res <- map_df(names(ts.exp)[-1],
                  Compare_means,ts.exp) %>% 
  arrange(gene)

DT::datatable(com_res)
## Only genes with a high pathway index are screened
com_res %>% 
  filter(pvalue<0.05) %>% 
  arrange(median.compare) %>% 
  filter(median.compare=="low < high") %>% 
  arrange(gene)->com_res_fin



com_res_fin %>% 
  arrange(gene) %>% 
  pull(gene) ->gene


library(org.Hs.eg.db)
cols <- c("GENENAME","ENTREZID")#Only the column containing the corresponding information in the package to be extracted by the gene with high pathway index is screened
select(org.Hs.eg.db,keys=gene, columns=cols, keytype="SYMBOL") %>% 
  drop_na() %>% 
  dplyr::slice(1:5) %>% 
  mutate(new=paste0(GENENAME," (",SYMBOL,")")) %>% 
  pull(new) %>% 
  paste(collapse = ", ")


library(patchwork)

creat_plot <- function(gene,dat=ts.exp){
dat %>% 
  select(gene,group) %>% 
  ggplot(aes_string("group",gene))+
  geom_violin(aes_string(fill = "group"), 
              draw_quantiles = c(0.25, 0.5, 0.75), 
              colour = "white", 
              alpha = 0.3, size = 0.8) +
  geom_boxplot(notch = TRUE, 
                width = 0.1, size = 0.8)+
  ggpubr::stat_compare_means(comparisons =list(c("high","low")) ,
                             size=5,
                             label = "p.signif")+
  theme_minimal(base_size = 18)+
  labs(x="Pathway index",y=gene)+
  theme(legend.position = "none")
    

    
  
}




if (T) {
  library(gghalves)
  library(ggpubr)
  halfbox <- function(gene,
                      dat,
                      group1color="#cba1d2",
                      group2color="#7067CF"){
    
    
 
    
    ## This is a very important step
    dat$groupnew <- jitter(as.numeric(factor(dat$group)), 
                           amount=.08)
    
    
    
    dat %>% 
      group_split(group) %>% 
      map_df(~data.frame(group=unique(.x$group),
                         mean=mean(.x[[gene]]), 
                         sd=sd(.x[[gene]])))->
      summary_df 
    
    dat$group = factor(as.numeric(factor(dat$group)))
    my_comparisons <- list(c("1","2"))
    
    ggplot(dat,aes_string(x="group",y=gene,fill="group",color="group")) +
      ## dot
      geom_point(aes(x=groupnew, colour=group),size =2) +
      ## Draw a line
      #geom_line(aes(x=groupnew,group=pairinfo), colour="black",alpha=0.1) +
      ## Draw the errorbar on the left
      # geom_errorbar(dat = summary_df[summary_df$group=="1",], 
      #               aes(x=group,y=mean, ymin=mean-sd, ymax=mean+sd), 
      #               position=position_nudge(x=.2),width = 0.05, size = 2,color="#cba1d2")+
      # geom_point(dat = summary_df[summary_df$group=="1",], aes(y=mean), position=position_nudge(x=.2), size = 10,color="#cba1d2")+
      # ## Draw the errorbar on the right
      # geom_errorbar(dat = summary_df[summary_df$group=="2",], 
      #               aes(x=group,y=mean, ymin=mean-sd, ymax=mean+sd), 
      #               position=position_nudge(x=-.2),width = 0.05, size = 2,color="#7067CF")+
    # geom_point(dat = summary_df[summary_df$group=="2",], aes(y=mean), position=position_nudge(x=-.2), size = 10,color="#7067CF")+
    # geom_line(dat = summary_df, aes(x = c(1.2,1.8), y = mean), color = 'black', size = 1) +
    # geom_point(dat = summary_df[summary_df$group=="1",], aes(y=mean), position=position_nudge(x=.2), size = 5,color="black")+
    # geom_point(dat = summary_df[summary_df$group=="2",], aes(y=mean), position=position_nudge(x=-.2), size = 5,color="black")+
    labs(x="Pathway index",y=gene) +
      theme_bw(base_size = 18)+
      theme(legend.position = "none")+
      ## Increment p-value
      stat_compare_means(comparisons = my_comparisons, 
                         paired=F,label = "p.signif",
                         size=5)+
      ## Left-hand box diagram
      geom_half_boxplot(data = dat[dat$group==1,],aes_string(x = "group", y=gene), position = position_nudge(x = -.25), side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2,fill = group1color,colour="black") +
      ## Right-hand box diagram
      geom_half_boxplot(data = dat[dat$group==2,],aes_string(x = "group", y=gene), position = position_nudge(x = .15), side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2,fill = group2color,colour="black") +
      ## The one on the left is the violin
      geom_half_violin(data = dat[dat$group==1,],aes_string(x = "group", y=gene), position = position_nudge(x = -.3), side = "l", fill = group1color,colour="black", width = .35) +
      ## On the right is a violin
      geom_half_violin(data = dat[dat$group==2,],aes_string(x = "group", y=gene), position = position_nudge(x = .3), side = "r", fill =  group2color,colour="black", width = .35) +
      
      ## Adjust the color of the fill
      scale_color_manual(values = c("1" = group1color, "2"=group2color)) +
      scale_fill_manual(values = c("1" = group1color, "2"=group2color))+
      ## Modify the X-axis label
      scale_x_discrete(labels=c("1" = "high", "2" = "low"))
    
  
}
}
ts.exp %>% 
  select(group,com_res_fin$gene)->
  dat.int

my.plot1 <- map(sort(names(dat.int)[-1]),creat_plot)

my.plot2 <- map(sort(names(dat.int)[-1]),
               halfbox,
               dat.int,
               group1color="#e54a34",
               group2color="#0e72ba")

cyt_plot <- halfbox("CYT",dat.int)
setwd("C:/Users/Cloudy/Desktop/论文写作/文章写作预分析/figure/表达谱箱图_修改")

ggsave(cyt_plot ,height = 5,
       width = 5,filename = "CYT_PLOT.pdf")


files = str_c(sort(names(dat.int)[-1]), "_exp.pdf") 


walk2(files,
  my.plot1,
  ggsave,
  height = 5,
  width = 5)
