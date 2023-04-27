## Analysis and visualization of immunoinfiltration results
library(tidyverse)
load("dat_input/TCGA NSCLC cohort Multivariate Cox analysis clean data.Rdata")
load("dat_input/NSCLC-exp-clin-data.Rdata")
load("dat_input/TCGA_NCSLC cohort after cibersort and ssGSEA.Rdata")
load("dat_output/Calculation result of TCGA cohort path index.Rdata")
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




cbind(LUAD_fpkmanno,LUSC_fpkmanno) %>% 
  t() %>% 
  as.data.frame() %>% 
  filter(substr(rownames(.),14,15)!="11") %>% 
  select(sort(pathway.up),sort(pathway.down))->
  gene.exp
  
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
  rownames_to_column("id") %>% 
  inner_join(index_clean %>%   #
               select(index_minus) %>% 
               rownames_to_column("id") 
             ) %>% 
  column_to_rownames("id") %>% 
  mutate(group=ifelse(index_minus>=0.108,"high","low")) %>% 
  select(-index_minus) %>% 
  select(group,everything())->
  gsva_index_nsclc


cell_com <- function(cell,dat){
    
   dat %>% 
   select(cell,group) %>% 
   Compare_means(cell,.)->res
   
     res
    
}


gsva_index_nsclc %>% 
  filter(substr(rownames(.),14,15)!="11") %>% 
  map_df(names(.)[-1],
         cell_com,
         dat=.) %>% 
  filter(pvalue<0.05&median.compare!="high=low")->
  cell.com.res


cell.com.res$median.compare <- factor(cell.com.res$median.compare,
                                       levels = c("low < high","high < low"))
cell.com.res <- arrange(  cell.com.res,median.compare)



gsva_index_nsclc %>% 
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
  mutate(group=ifelse(index_minus>=0.108,"high","low")) %>% 
  select(-index_minus) %>% 
  select(group,everything())->
  ciber_index_nsclc


ciber_index_nsclc %>% 
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

ciber_index_nsclc  %>% 
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



c
d
x11()



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
  mutate(group=ifelse(index_minus>=0.108,"high","low")) %>% 
  select(group,everything()) %>% 
  select("GZMA","PRF1",any_of(gene_extend),"index_minus")->
  dat_int_gene

dat_int_gene %>% 
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


# Statistical result
com_res <- map_df(names(ts.exp)[-1],
                  Compare_means,ts.exp) %>% 
  arrange(gene)

## Only genes with a high pathway index are screened
com_res %>% 
  filter(pvalue<0.05) %>% 
  arrange(median.compare) %>% 
  filter(median.compare=="low < high") %>% 
  arrange(gene)->com_res_fin

com_res_fin %>% 
  arrange(gene) %>% 
  pull(gene) ->gene




# ploting -----------------------------------------------------------------




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



if (F) {
  
  
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
ggsave(cyt_plot ,height = 5,
       width = 5,filename = "CYT_PLOT.pdf")


files = str_c(sort(names(dat.int)[-1]), "_exp.pdf") 


walk2(files,
  my.plot1,
  ggsave,
  height = 5,
  width = 5)
