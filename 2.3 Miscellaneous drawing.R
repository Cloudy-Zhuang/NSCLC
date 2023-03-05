library(tidyverse)
library(ggpubr)
##Data from the clinical matrix were incorporated
load("dat_input/NSCLC-exp-clin-data.Rdata")
load("dat_output/TCGA队列通路指数计算结果.Rdata")


NSCLC_fpkmanno <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)
NSCLC_clin_tab <- rbind(LUAD_clin_tab,LUSC_clin_tab)

Creat_dat.os <- function(fpkmanno,clin_tab,mygene,
                         survival){
  
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
  
}

index_dat <- as.data.frame(t(index_clean))

Creat_dat.os(fpkmanno =index_dat,
             mygene = "index_minus",
             NSCLC_clin_tab) ->
  dat.os.multicox
dat.os.multicox[is.na(dat.os.multicox)] <- "Unknown"
save(dat.os.multicox,file="dat_output/多因素cox分析清洁数据.Rdata")

#Tide score
tide_low <- read.csv("dat_input/index_low TIDE.csv") %>% 
  mutate(group="TIDE_low")
tide_high <- read.csv("dat_input/index_high TIDE.csv") %>% 
  mutate(group="TIDE_high")

tide_high <-read.csv("dat_input/不校正的高指数组TIDE评分.csv") %>% 
  mutate(group="high")
tide_low <-read.csv("dat_input/不校正低指数组TIDE评分.csv") %>% 
  mutate(group="low")

load("dat_output/多因素cox分析清洁数据.Rdata")


rbind(tide_low,tide_high) %>% 
  filter(substr(Patient,14,15)!="11") %>% 
  #mutate(Patient=substr(Patient,1,12)) %>% 
  #filter(Patient%in%rownames(dat.os.multicox)) %>% 
  ggplot(aes(group,TIDE,color=group))+
  geom_boxplot()+
  geom_jitter(alpha = 0.5)+
  stat_compare_means(size=5,
                     label.x = 1.3,
                     show.legend = F,
                     comparisons = list(c("high","low")),
                     label = "p.signif")+
  ggsci::scale_color_jco()+
  labs(x="Pathway index",
       y="TIDE Score")+
  ggsci::scale_color_nejm()+
  theme_bw(base_size = 18)+
  theme(legend.position = "none")->p
x11()
ggsave(p,filename = "TIDE评分修正图.pdf",
       height = 5,
       width = 5)


rbind(tide_low,tide_high) %>% 
  mutate(Patient=substr(Patient,1,12)) %>%
  distinct(Patient,.keep_all = T) %>% 
  inner_join(dat.os.multicox %>% 
               rownames_to_column("Patient") %>% 
               select(Patient,index_minus))->
  patient_index
x11()
cor.test(patient_index$TIDE,patient_index$index_minus)

#Panpie
##Pick five factors
dat.os.multicox


# age
# Survival outcome
# Treatment outcome
#Clinical staging
# Receive radiation therapy
dat<- dat.os.multicox %>% 
           mutate(Risk=ifelse(index_minus>=0.108,"High","Low")) %>% 
           select(Risk,
                  age,
                  surv_status,
                  Stage,
                  radiation_therapy,
                  Treatment_outcome
                  ) 
dat$radiation_therapy[dat$radiation_therapy==""] <- "Unknown"
openxlsx::write.xlsx(dat,"饼图示例数据.xlsx")
x11()
if (T) {
  # 按Risk分成High和Low，Calculate each column value。
  gname <- "Risk"
  vname <- setdiff(colnames(dat), gname)
  pie.high <- pie.low <- list()
  fisher.p <- c()
  for (i in vname) {
    
    tmp <- table(dat[,gname], dat[,i])
    p <- format(fisher.test(tmp)$p.value,digits = 2)
    names(p) <- i
    fisher.p <- c(fisher.p, p)
    
    pie.dat <- 
      tmp %>% as.data.frame() %>%
      group_by(Var1) %>%
      mutate(Pct = Freq/sum(Freq)) %>% 
      as.data.frame()
    
    # The two rows in the table correspond to the two categories of Risk: Risk high and Risk low
    pie.high[[i]] <- pie.dat[which(pie.dat$Var1 == "High"),]
    pie.low[[i]] <- pie.dat[which(pie.dat$Var1 == "Low"),]
  }
  
  
  # Set color
  black  <- "#1E1E1B"
  blue   <- "#3C4E98"
  yellow <- "#E4DB36"
  orange <- "#E19143"
  green  <- "#57A12B"
  cherry <- "#8D3A86"
  
  # Create color
  status.col <- c("grey80",black)
  stage.col <- alpha(blue, c(0.4, 0.6, 0.8, 1))
  M.col <- c(yellow, orange)
  #N.col <- alpha(green, c(0.5, 0.7, 1))
  N.col <- c("#fb6501","#6699cc")
  #T.col <- alpha(cherry, c(0.4, 0.6, 0.8, 1))
  T.col <- c("#19caad", "#8cc7b5","#a0eee1","#beedc7")
 
   #The core base plot can be drawn one piece at a time. Of course, pie chart can also be extracted for later AI or PPT splicing
  #pdf("pieTable.pdf",width = 7, height = 5)
  showLayout <- T# By default, the layout structure is not displayed on the first page of the final pdf, but it is recommended to change to TRUE when the first drawing, easy to understand
  
  # Set the layout of the screen, the same number represents the same block, the more number represents the larger the area of the block (a total of 25 areas)
  
  layout(matrix(c( 1, 1, 1,  2, 2, 2,  3, 3, 3,  4, 4, 4,  5, 5, 5,  6, 6, 6,
                   7, 7, 7,  8, 8, 8,  9, 9, 9, 10,10,10, 11,11,11, 12,12,12,
                   7, 7, 7,  8, 8, 8,  9, 9, 9, 10,10,10, 11,11,11, 12,12,12,
                   13,13,13, 14,14,14, 15,15,15, 16,16,16, 17,17,17, 18,18,18,
                   13,13,13, 14,14,14, 15,15,15, 16,16,16, 17,17,17, 18,18,18,
                   19,19,19, 20,20,20, 21,21,21, 22,22,22, 23,23,23, 24,24,24,
                   25,25,25, 25,25,25, 25,25,25, 25,25,25, 25,25,25, 25,25,25,),
                byrow = T,nrow = 7))
  
  if(showLayout) {
    layout.show(n = 25) # Visual display canvas distribution
  }
  
  #-------------------------#
  # Canvas area 1-6: Draw the header #
  #-------------------------#
  
  par(bty="n", mgp = c(0,0,0), mar = c(0,0,0,0), lwd = 2) # Base parameter, the distance of each boundary is 0
  plot(1,1,
       xlab = "",xaxt = "n", # The X-axis is not displayed
       ylab = "",yaxt = "n") # The Y-axis is not displayed
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # Background blackening
  text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")To get the absolute position of the canvas
       (par("usr")[3]+par("usr")[4])/2,
       "Pathway\nindex",cex = 2, col = "white") # Display title
  
  plot(1,1,
       xlab = "",xaxt = "n", # The X-axis is not displayed
       ylab = "",yaxt = "n") # The Y-axis is not displayed
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # Background blackening
  text((par("usr")[1]+par("usr")[2])/2, # use par("usr")To get the absolute position of the canvas
       (par("usr")[3]+par("usr")[4])/2,
       "Age",cex = 2, col = "white") # Display title
  
  plot(1,1,
       xlab = "",xaxt = "n", # The X-axis is not displayed
       ylab = "",yaxt = "n") # The Y-axis is not displayed
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # Background blackening
  text((par("usr")[1]+par("usr")[2])/2, # usepar("usr")To get the absolute position of the canvas
       (par("usr")[3]+par("usr")[4])/2,
       "Survival\nStatus",cex = 2, col = "white") # Display title
  
  plot(1,1,
       xlab = "",xaxt = "n", # The X-axis is not displayed
       ylab = "",yaxt = "n") # The Y-axis is not displayed
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # Background blackening
  text((par("usr")[1]+par("usr")[2])/2, # usepar("usr")To get the absolute position of the canvas
       (par("usr")[3]+par("usr")[4])/2,
       "Stage",cex = 2, col = "white") # Display title
  
  plot(1,1,
       xlab = "",xaxt = "n", # The X-axis is not displayed
       ylab = "",yaxt = "n") # The Y-axis is not displayed
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # Background blackening
  text((par("usr")[1]+par("usr")[2])/2, # usepar("usr")To get the absolute position of the canvas
       (par("usr")[3]+par("usr")[4])/2,
       "Radition\nTheraphy",cex = 2, col = "white") # Display title
  
  plot(1,1,
       xlab = "",xaxt = "n", # The X-axis is not displayed
       ylab = "",yaxt = "n") # The Y-axis is not displayed
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # Background blackening
  text((par("usr")[1]+par("usr")[2])/2, # use par("usr")To get the absolute position of the canvas
       (par("usr")[3]+par("usr")[4])/2,
       "Treatment\nOutcome",cex = 2, col = "white") # Display title
  
  #--------------------------------------#
  # Canvas area 7-12: Draw head and fan charts for High group #
  #--------------------------------------#
  
  plot(1,1,
       xlab = "",xaxt = "n", # The X-axis is not displayed
       ylab = "",yaxt = "n") # The Y-axis is not displayed
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # Background blackening
  text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")To get the absolute position of the canvas
       (par("usr")[3]+par("usr")[4])/2,
       "High\n(n = 858)",cex = 2, col = "white") # Display title
  
  # High group
  pie(pie.high$age$Pct, 
      col = status.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  
  pie(pie.high$surv_status$Pct, 
      col = stage.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  
  pie(pie.high$Stage$Pct, 
      col = M.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  
  pie(pie.high$radiation_therapy$Pct, 
      col = N.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  
  pie(pie.high$Treatment_outcome$Pct, 
      col = T.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  abline(v = par("usr")[2], col = "black") # Seal the black line on the right side
  
  #--------------------------------------#
  # Canvas area 13-18: Draw the Low group head and fan chart #
  #--------------------------------------#
  
  plot(1,1,
       xlab = "",xaxt = "n", # The X-axis is not displayed
       ylab = "",yaxt = "n") # The Y-axis is not displayed
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # Background blackening
  text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")To get the absolute position of the canvas
       (par("usr")[3]+par("usr")[4])/2,
       "Low\n(n = 137)",cex = 2, col = "white") # Display title

  
  # Low group
  pie(pie.low$age$Pct, 
      col = status.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  
  pie(pie.low$surv_status$Pct, 
      col = stage.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  
  pie(pie.low$Stage$Pct, 
      col = M.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  
  pie(pie.low$radiation_therapy$Pct, 
      col = N.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  
  pie(pie.low$Treatment_outcome$Pct, 
      col = T.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  abline(v = par("usr")[2], col = "black") # Seal the black line on the right side
  
  #--------------------------------#
  # Canvas Area 19-24: Draw empty headers and p values
 #
  #--------------------------------#
  
  plot(1,1,
       xlab = "",xaxt = "n", # The X-axis is not displayed
       ylab = "",yaxt = "n") # The Y-axis is not displayed
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # Background blackening
  
  plot(1,1,col = "white",
       xlab = "",xaxt = "n", # The X-axis is not displayed
       ylab = "",yaxt = "n") # The Y-axis is not displayed
  text((par("usr")[1]+par("usr")[2])/2, # Use par("usr") to get the absolute location of the canvas
       (par("usr")[3]+par("usr")[4])/2,
       paste0("p = ",fisher.p["age"]),cex = 1.5, col = "black") # Display title
  abline(h = par("usr")[3], col = "black") # Seal the bottom with black lines
  
  plot(1,1,col = "white",
       xlab = "",xaxt = "n", # The X-axis is not displayed
       ylab = "",yaxt = "n") # The Y-axis is not displayed
  text((par("usr")[1]+par("usr")[2])/2, # Use par("usr") to get the absolute location of the canvas
       (par("usr")[3]+par("usr")[4])/2,
       paste0("p = ",fisher.p["surv_status"]),cex = 1.5, col = "black") # Display title
  abline(h = par("usr")[3], col = "black")
  
  plot(1,1,col = "white",
       xlab = "",xaxt = "n",# The X-axis is not displayed
       ylab = "",yaxt = "n") # The Y-axis is not displayed
  text((par("usr")[1]+par("usr")[2])/2, # Use par("usr") to get the absolute location of the canvas
       (par("usr")[3]+par("usr")[4])/2,
       paste0("p = ",fisher.p["Stage"]),cex = 1.5, col = "black") # Display title
  abline(h = par("usr")[3], col = "black")
  
  plot(1,1,col = "white",
       xlab = "",xaxt = "n", # The X-axis is not displayed
       ylab = "",yaxt = "n") # The Y-axis is not displayed
  text((par("usr")[1]+par("usr")[2])/2, # Use par("usr") to get the absolute location of the canvas
       (par("usr")[3]+par("usr")[4])/2,
       paste0("p = ",fisher.p["radiation_therapy"]),cex = 1.5, col = "black") # Display title
  abline(h = par("usr")[3], col = "black")
  
  plot(1,1,col = "white",
       xlab = "",xaxt = "n", # The X-axis is not displayed
       ylab = "",yaxt = "n") # The Y-axis is not displayed
  text((par("usr")[1]+par("usr")[2])/2, # Use par("usr") to get the absolute location of the canvas
       (par("usr")[3]+par("usr")[4])/2,
       paste0("p = ",fisher.p["Treatment_outcome"]),cex = 1.5, col = "black") # Display title
  abline(h = par("usr")[3], col = "black") # Seal the bottom with black lines
  abline(v = par("usr")[2], col = "black") # Seal the black line on the right side
  
  #----------------------#
  # Canvas Area 25: Draw the legend #
  #----------------------#
  
  plot(0,0,col = "white",
       xlab = "",xaxt = "n", # The X-axis is not displayed
       ylab = "",yaxt = "n") # The Y-axis is not displayed
  
  legend("topleft",
         legend = c("<65","≥65",
                    "Alive","Dead",
                    "I-II","III-IV",
                    "NO","YES",
                    "PD","SD","PR","CR"),
         fill = c(status.col,
                  stage.col[1:2],
                  M.col,
                  N.col,
                  T.col),
         border = NA, # Legend colors have no border
         bty = "n", # Legend has no border
         cex = 1.8,
         box.lwd = 3,
         x.intersp = 0.05,
         y.intersp = 1,
         text.width = 0.075, # The interval of the legend
         horiz = T) # The legend is placed horizontally
  
  
}




x11()



#重要基因23个基因的表达情况

load("dat_output/Mtor通路研究基因.Rdata")
load("dat_input/GTex_fpkm_anno_pheno.Rdata")
load("dat_output/通路基因生存分析及表达谱分析结果.Rdata")

NSCLC_clin <- rbind(LUAD_clin_tab,LUSC_clin_tab)
NSCLC_fpkm <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)


# Identification and matching of mtor gene ------------------------------------------------------------------

str_replace(mtor_hubgene,"GRK2","ADRBK1") %>% n_distinct()

interest_gene <- sort(unique(c(mtor_hubgene,motor_gene)))
interest_gene <- sort(replace(interest_gene,interest_gene=="GRK2","ADRBK1"))

# Pathway-related genes in the literature
gene_publish <- c("PIK3R1", "PIK3CA", "PIK3C2A", "PIK3C3", "AKT1", "PIK3R2", "PIK3CB", 
                  "PIK3C2B", "PIK3R4", "AKT2", "PIK3R3", "PIK3CD", "PIK3C2G", "PTEN", 
                  "AKT3", "TSC2", "MTOR","RHEB",
                  "PDK1","INPP4B","PHLPP1","STK11",
                  "TSC1")
gene_publish <- sort(gene_publish)
gene_publish %>% n_distinct()
gene.all <- sort(unique(c(interest_gene,gene_publish)))

## Exploration in the transcriptome
NSCLC_fpkm <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)
ts <- NSCLC_fpkm[1:20]
NSCLC_fpkm %>% 
  t() %>% 
  as.data.frame() %>% 
  select(gene.all) %>% 
  mutate(group=if_else(substr(names(NSCLC_fpkm),14,15)=="11",
                       "TCGA-normal","TCGA-cancer")) %>% 
  pivot_longer(-group,names_to = "gene",values_to = "exp")->
  NSCLC_fpkm_gene.all


# Normal lung tissue data
GTex_fpkm_anno_pheno %>% 
  mutate(across(where(is.numeric),~log2(((.x^2)-0.001)+1))) %>% 
  filter(tissue=="Lung") %>% 
  mutate(group="GTEx-normal") %>% 
  select(group,gene.all)->
  GTex_fpkm_lung

GTex_fpkm_lung %>% 
  pivot_longer(-group,names_to = "gene",values_to = "exp")->
  GTex_fpkm_gene.all


gene_exp.all <- rbind(NSCLC_fpkm_gene.all,GTex_fpkm_gene.all)

gene.important_sel <- setdiff(gene.important,"RPS6KA3") %>% 
  sort()

gene_exp.all %>% 
  filter(gene%in%gene.important_sel) ->
  gene.imp_sel

setwd("C:/Users/Cloudy/Desktop/论文写作/文章写作预分析/figure/表达谱箱图")


creat.box <- function(mygene){
gene.imp_sel %>% 
  filter(gene==mygene) %>% 
  mutate(group=factor(group,levels = unique(gene.imp_sel$group))) %>% 
  ggplot(aes(group,exp,fill=group))+
  geom_boxplot(width = 0.5)+
  labs(x=mygene,y=paste0(mygene," Expression"))+
  ggprism::theme_prism()+
  theme(legend.position = "none")+
  scale_fill_manual(values = c( "#FC4E07", "#E7B800","#00AFBB"))+
  stat_compare_means(comparisons = combn(unique(gene.imp_sel$group),2,simplify = F),
                     label = "p.signif",size=5)->box 

  ggsave(box ,filename = paste0(mygene,"_boxplot.pdf"),
         height = 5,width = 5)
}

walk(gene.important_sel,creat.box)


## Single gene survival analysis forest map




##17 pathway gene survival analysis Rose Map


##
index_clean %>% 
  mutate(group=ifelse(substr(rownames(.),14,15)!="11","cancer","para-cancer")) %>% 
  ggplot(aes(group,index_minus,fill=group))+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.75)+
  geom_boxplot(width=0.2,
               size=0.75,outlier.colour = NA )+
  ggpubr::stat_compare_means(show.legend = F,size=5,
                             label.x = 1.2)+
  theme_bw(base_size = 18)+
  theme(legend.position = "none")+
  labs(y="Pathway index") ->p
  ggsave(p,filename = "pi_incancer.pdf",width = 5,height = 5)


  
  

# ROC curve -------------------------------------------------------------------
  index_clean %>% 
    mutate(group=ifelse(substr(rownames(.),14,15)!="11","cancer","para-cancer")) %>% 
    select(index_minus,group)->ts
  
  roc_plot("index_minus",data=ts,outcome = "group")
  