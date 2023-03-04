library(tidyverse)
library(gtsummary)


# load data ---------------------------------------------------------------
load("dat_input/NSCLC-exp-clin-data.Rdata")
mut.dat <- data.table::fread("dat_input/mc3.v0.2.8.PUBLIC.xena.gz",data.table = F)


NSCLC_clin <- rbind(LUAD_clin_tab,LUSC_clin_tab)
NSCLC_fpkm <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)
mut.dat %>% 
  select(sample,effect,gene) %>%
  filter(sample%in%substr(names(NSCLC_fpkm),1,15)) %>% 
  group_split(sample,gene) %>% 
  map_df(~data.frame(sample=unique(.x$sample),
                     gene=unique(.x$gene),
                     effect=paste(.x$effect,collapse = "/")))->
  mut.dat_new

mut.dat_new$gene[grep("MET",mut.dat_new$gene)] %>% unique()



NSCLC_mut.dat <- mut.dat_new

NSCLC_mut.dat %>% 
  filter(gene%in%sort(c("EGFR","ALK","NTRK1","ROS1","RET","NTRK2","NTRK3","ERBB2",
                        "BRAF","RET","KRAS","MAP2K1","FLT1","MET","MST1R"))) %>% 
  group_by(gene) %>% 
  mutate(index=row_number()) %>% 
  pivot_wider(names_from = gene,
              values_from = effect ) %>% 
  dplyr::select(sample,any_of(sort(c("EGFR","ALK","NTRK1","ROS1","RET","NTRK2",
                              "NTRK3","BRAF","RET","KRAS","MAP2K1",
                              "FLT1"))))->
  NSCLC_sen.mut

load("dat_output/TCGA队列通路指数计算结果.Rdata")
NSCLC_sen.mut %>% 
  mutate(across(-1,function(x){x=ifelse(is.na(x),"No","Yes")})) %>%
  pivot_longer(cols = -1,names_to = "gene_mut",values_to = "mut_status") %>% 
  inner_join(index_clean %>% 
               dplyr::select(index_minus) %>% 
               rownames_to_column("sample") %>% 
               mutate(sample=substr(sample,1,15)) %>% 
               filter(substr(sample,14,15)!="11")
  ) %>% 
  distinct() %>% 
  mutate(mut_status=factor(mut_status,levels = c("Yes","No")))%>% 
  ggplot(aes(x = gene_mut, y = index_minus))+
  geom_boxplot(aes(fill = mut_status),
               position = position_dodge(1),width=.3,outlier.shape = NA)+
  # geom_violin(aes(colour = mut_status),position = position_dodge(1),scale = "width",fill=NA)+
  theme_bw(base_size = 18)+
  labs(fill="Mutation status",
       color="Mutation status",
       x="Targeted therapy-related genes",
       y="Pathway index")+
  ggsci::scale_color_npg()+
  ggsci::scale_fill_npg()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                   vjust = 1, colour = "black"),
        legend.position = "top")+
  ggpubr::stat_compare_means(aes(group=mut_status), 
                             label = "p.signif")+
  geom_hline(yintercept = 0.108,linetype=2,col="red",size=0.8,alpha=0.5)+
  scale_y_continuous(breaks = c(-0.5,0,0.108,0.5))+
  x11()




#save(mut.dat_new,NSCLC_sen.mut,file = "dat_output/TCGA突变数据.Rdata")

if (F) {
  
  NSCLC_clin$ECOG
  
  NSCLC_clin %>% 
    filter(rownames(.)%in%substr(names(NSCLC_fpkm),1,12)) %>% 
    mutate_all(as.character) %>% 
    na_if("") %>% 
    mutate_all(replace_na,"Unknown") %>% 
    mutate(post_bronchodilator_fev1_fvc_percent=as.numeric(post_bronchodilator_fev1_fvc_percent),
           pack_years_smoked.exposures=as.numeric(pack_years_smoked.exposures),
           OS.time=as.numeric(OS.time)/365,
           post_bronchodilator_fev1_percent=as.numeric(post_bronchodilator_fev1_percent)) %>% 
    tbl_summary(
      statistic=list(all_continuous() ~ "{mean}±{sd} ({p25}, {median},{p75},)"),
      sort = list(everything() ~ "alphanumeric" ),
    )
}




  
if (F) {
  mygene <- c("EGFR","ALK","NTRK1","ROS1","RET","NTRK2",
              "NTRK3","BRAF","RET","KRAS","MAP2K1",
              "FLT1") %>% sort()
  
  
  
  
  library(org.Hs.eg.db)
  mygene %>% 
    paste(collapse = ", ")
  
  ##For example, you want to get its abbreviation and ENTREZID through the full name of the gene. The gene name here is keytype
  ensids <- c("tumor protein p53")#Full name of gene
  cols <- "GENENAME"#The column of the corresponding information contained in the package you want to extract
  select(org.Hs.eg.db, keys=mygene, columns=cols, keytype="SYMBOL") %>% 
    as.data.frame() %>% 
    mutate(gene=paste0(GENENAME," (",SYMBOL,")")) %>% 
    pull(gene) %>% 
    paste(collapse = ", ")
}




NSCLC_sen.mut %>% 
  mutate(across(-1,function(x){x=ifelse(is.na(x),"No","Yes")})) %>%
  inner_join(index_clean %>% 
               dplyr::select(index_minus) %>% 
               rownames_to_column("sample") %>% 
               mutate(sample=substr(sample,1,15)) %>% 
               filter(substr(sample,14,15)!="11")) %>% 
  mutate(Risk=ifelse(index_minus>0.108,"High","Low")) %>% 
  dplyr::select(-index_minus,-sample) %>%
  dplyr::select(Risk,everything())->ts
  
library(gtsummary)
ts %>% 
  tbl_summary(by =Risk ) %>% 
  add_p(all_categorical() ~ "fisher.test")




if (T) {
  
  dat <- ts %>% 
    select(-"MAP2K1") #Delete another meaningless one and keep 10
  dat$Risk %>% table()
  dat_onetosix=dplyr::select(dat,1,7:11) #5 genes are input once
  x11()
  if (T) {
    # Divide the risk into High and Low, and calculate the values of each column.
    gname <- "Risk"
    vname <- setdiff(colnames(dat_onetosix), gname)
    pie.high <- pie.low <- list()
    fisher.p <- c()
    for (i in vname) {
      dat0 <-dat_onetosix %>%  
        dplyr::select(gname,i) %>% 
        drop_na()
      
      tmp <- table(dat0[,gname,drop=T], dat0[,i,drop=T])
      
      
      p <- format(fisher.test(tmp)$p.value,digits = 2)
      
      
      names(p) <- vname[which(vname==i)]
      fisher.p <- c(fisher.p, p)
      
      pie.dat <- 
        tmp %>% as.data.frame() %>%
        group_by(Var1) %>%
        mutate(Pct = Freq/sum(Freq)) %>% 
        as.data.frame()
      
      # The two rows in the table correspond to two types of risk: risk high and risk low
      pie.high[[i]] <- pie.dat[which(pie.dat$Var1 == "High"),]
      pie.low[[i]] <- pie.dat[which(pie.dat$Var1 == "Low"),]
    }
    
    
    # Set Color
    black  <- "#1E1E1B"
    blue   <- "#3C4E98"
    yellow <- "#E4DB36"
    orange <- "#E19143"
    green  <- "#57A12B"
    cherry <- "#8D3A86"
    
    # Create Color
    status.col <- c("grey80",black)
    stage.col <- alpha(blue, c(0.4, 0.6, 0.8, 1))
    M.col <- c(yellow, orange)
    #N.col <- alpha(green, c(0.5, 0.7, 1))
    N.col <- c("#fb6501","#6699cc")
    #T.col <- alpha(cherry, c(0.4, 0.6, 0.8, 1))
    T.col <- c("#19caad", "#8cc7b5","#a0eee1","#beedc7")
    
    # The hard core base plot is painted one by one. Of course, it is also convenient to extract the pie chart from it for later AI or PPT splicing
    #pdf("pieTable.pdf",width = 7, height = 5)
    showLayout <- T# By default, the layout structure is not displayed on the front page of the final pdf. However, it is recommended to change it to TRUE when drawing for the first time for easy understanding
    
    # Set the screen layout. The same number represents the same block. The more the number, the larger the area of the block (a total of 25 areas)
    
    layout(matrix(c( 1, 1, 1,  2, 2, 2,  3, 3, 3,  4, 4, 4,  5, 5, 5,  6, 6, 6,
                     7, 7, 7,  8, 8, 8,  9, 9, 9, 10,10,10, 11,11,11, 12,12,12,
                     7, 7, 7,  8, 8, 8,  9, 9, 9, 10,10,10, 11,11,11, 12,12,12,
                     13,13,13, 14,14,14, 15,15,15, 16,16,16, 17,17,17, 18,18,18,
                     13,13,13, 14,14,14, 15,15,15, 16,16,16, 17,17,17, 18,18,18,
                     19,19,19, 20,20,20, 21,21,21, 22,22,22, 23,23,23, 24,24,24,
                     25,25,25, 25,25,25, 25,25,25, 25,25,25, 25,25,25, 25,25,25),
                  byrow = T,nrow = 7))
    
    if(showLayout) {
      layout.show(n = 25) # Visual display of canvas distribution
    }
    
    #-------------------------------#
    # Canvas area 1-6: drawing head #
    #-------------------------------#
    
    par(bty="n", mgp = c(0,0,0), mar = c(0,0,0,0), lwd = 2) # Basic parameters, each boundary distance is 0
    plot(1,1,
         xlab = "",xaxt = "n", # Do not display the X-axis
         ylab = "",yaxt = "n") # Do not display the Y-axis
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # Black background
    text((par("usr")[1]+par("usr")[2])/2, # Use par ("usr") to obtain the absolute position of the canvas
         (par("usr")[3]+par("usr")[4])/2,
         "Pathway\nindex",cex = 2, col = "white") # Display Chart Title
    
    plot(1,1,
         xlab = "",xaxt = "n", # Do not display the X-axis
         ylab = "",yaxt = "n") # Do not display the Y-axis
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # Black background
    text((par("usr")[1]+par("usr")[2])/2, # Use par ("usr") to obtain the absolute position of the canvas
         (par("usr")[3]+par("usr")[4])/2,
         vname[1],cex = 2, col = "white") # Display Chart Title
    
    plot(1,1,
         xlab = "",xaxt = "n", # Do not display the X-axis
         ylab = "",yaxt = "n") # Do not display the Y-axis
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # Black background
    text((par("usr")[1]+par("usr")[2])/2, # Use par ("usr") to obtain the absolute position of the canvas
         (par("usr")[3]+par("usr")[4])/2,
         vname[2],cex = 2, col = "white") # Display Chart Title
    
    plot(1,1,
         xlab = "",xaxt = "n", # Do not display the X-axis
         ylab = "",yaxt = "n") # Do not display the Y-axis
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # Black background
    text((par("usr")[1]+par("usr")[2])/2, # Use par ("usr") to obtain the absolute position of the canvas
         (par("usr")[3]+par("usr")[4])/2,
         vname[3],cex = 2, col = "white") # Display Chart Title
    
    plot(1,1,
         xlab = "",xaxt = "n", # Do not display the X-axis
         ylab = "",yaxt = "n") # Do not display the Y-axis
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # Black background
    text((par("usr")[1]+par("usr")[2])/2, # Use par ("usr") to obtain the absolute position of the canvas
         (par("usr")[3]+par("usr")[4])/2,
         vname[4],cex = 2, col = "white") # Display Chart Title
    
    plot(1,1,
         xlab = "",xaxt = "n", # Do not display the X-axis
         ylab = "",yaxt = "n") # Do not display the Y-axis
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # Black background
    text((par("usr")[1]+par("usr")[2])/2, # Use par ("usr") to obtain the absolute position of the canvas
         (par("usr")[3]+par("usr")[4])/2,
         vname[5],cex = 2, col = "white") # Display Chart Title
    
    #---------------------------------------------------#
    # Canvas area 7-12: Draw high group head and sector #
    #---------------------------------------------------#
    
    plot(1,1,
         xlab = "",xaxt = "n", # Do not display the X-axis
         ylab = "",yaxt = "n") # Do not display the Y-axis
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # Black background
    text((par("usr")[1]+par("usr")[2])/2, # Use par ("usr") to obtain the absolute position of the canvas
         (par("usr")[3]+par("usr")[4])/2,
         "High\n(n = 726)",cex = 2, col = "white") # Display Chart Title
    
    # High group
    pie(pie.high[[1]]$Pct, 
        col = c("#16a085","#e74c3c"), 
        border = "white",  
        radius = 1, 
        labels = NA,
        init.angle = 90)
    symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
    
    pie(pie.high[[2]]$Pct, 
        col =c("#16a085","#e74c3c"), 
        border = "white",  
        radius = 1, 
        labels = NA,
        init.angle = 90)
    symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
    
    pie(pie.high[[3]]$Pct, 
        col = c("#16a085","#e74c3c"), 
        border = "white",  
        radius = 1, 
        labels = NA,
        init.angle = 90)
    symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
    
    pie(pie.high[[4]]$Pct, 
        col = c("#16a085","#e74c3c"), 
        border = "white",  
        radius = 1, 
        labels = NA,
        init.angle = 90)
    symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
    
    pie(pie.high[[5]]$Pct, 
        col = c("#16a085","#e74c3c"), 
        border = "white",  
        radius = 1, 
        labels = NA,
        init.angle = 90)
    symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
    abline(v = par("usr")[2], col = "black") # Seal the black line on the right side
    
    #---------------------------------------------------#
    # Canvas area 13-18: Draw low group head and sector #
    #---------------------------------------------------#
    
    plot(1,1,
         xlab = "",xaxt = "n", # Do not display the X-axis
         ylab = "",yaxt = "n") # Do not display the Y-axis
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # Black background
    text((par("usr")[1]+par("usr")[2])/2, # Use par ("usr") to obtain the absolute position of the canvas
         (par("usr")[3]+par("usr")[4])/2,
         "Low\n(n = 104)",cex = 2, col = "white") # Display Chart Title
    
    # Low group
    pie(pie.low[[1]]$Pct, 
        col = c("#16a085","#e74c3c"), 
        border = "white",  
        radius = 1, 
        labels = NA,
        init.angle = 90)
    symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
    
    pie(pie.low[[2]]$Pct, 
        col = c("#16a085","#e74c3c"), 
        border = "white",  
        radius = 1, 
        labels = NA,
        init.angle = 90)
    symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
    
    pie(pie.low[[3]]$Pct, 
        col = c("#16a085","#e74c3c"), 
        border = "white",  
        radius = 1, 
        labels = NA,
        init.angle = 90)
    symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
    
    pie(pie.low[[4]]$Pct, 
        col = c("#16a085","#e74c3c"), 
        border = "white",  
        radius = 1, 
        labels = NA,
        init.angle = 90)
    symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
    
    pie(pie.low[[5]]$Pct, 
        col = c("#16a085","#e74c3c"), 
        border = "white",  
        radius = 1, 
        labels = NA,
        init.angle = 90)
    symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
    abline(v = par("usr")[2], col = "black") # Seal the black line on the right side
    
    #------------------------------------------------#
    # Canvas area 19-24: Draw empty head and p value #
    #------------------------------------------------#
    
    plot(1,1,
         xlab = "",xaxt = "n", # Do not display the X-axis
         ylab = "",yaxt = "n") # Do not display the Y-axis
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # Black background
    
    plot(1,1,col = "white",
         xlab = "",xaxt = "n", # Do not display the X-axis
         ylab = "",yaxt = "n") # Do not display the Y-axis
    text((par("usr")[1]+par("usr")[2])/2, # Use par ("usr") to obtain the absolute position of the canvas
         (par("usr")[3]+par("usr")[4])/2,
         paste0("p = ",fisher.p[vname[1]]),cex = 1.5, col = "black") # Display Chart Title
    abline(h = par("usr")[3], col = "black") # Seal the black line at the bottom
    
    plot(1,1,col = "white",
         xlab = "",xaxt = "n", # Do not display the X-axis
         ylab = "",yaxt = "n") # Do not display the Y-axis
    text((par("usr")[1]+par("usr")[2])/2, # Use par ("usr") to obtain the absolute position of the canvas
         (par("usr")[3]+par("usr")[4])/2,
         paste0("p = ",fisher.p[vname[2]]),cex = 1.5, col = "black") # Display Chart Title
    abline(h = par("usr")[3], col = "black")
    
    plot(1,1,col = "white",
         xlab = "",xaxt = "n",# Do not display the X-axis
         ylab = "",yaxt = "n") # Do not display the Y-axis
    text((par("usr")[1]+par("usr")[2])/2, # Use par ("usr") to obtain the absolute position of the canvas
         (par("usr")[3]+par("usr")[4])/2,
         paste0("p = ",fisher.p[vname[3]]),cex = 1.5, col = "black") # Display Chart Title
    abline(h = par("usr")[3], col = "black")
    
    plot(1,1,col = "white",
         xlab = "",xaxt = "n", # Do not display the X-axis
         ylab = "",yaxt = "n") # Do not display the Y-axis
    text((par("usr")[1]+par("usr")[2])/2, # Use par ("usr") to obtain the absolute position of the canvas
         (par("usr")[3]+par("usr")[4])/2,
         paste0("p = ",fisher.p[vname[4]]),cex = 1.5, col = "black") # Display Chart Title
    abline(h = par("usr")[3], col = "black")
    
    plot(1,1,col = "white",
         xlab = "",xaxt = "n", # Do not display the X-axis
         ylab = "",yaxt = "n") # Do not display the Y-axis
    text((par("usr")[1]+par("usr")[2])/2, # Use par ("usr") to obtain the absolute position of the canvas
         (par("usr")[3]+par("usr")[4])/2,
         paste0("p = ",fisher.p[vname[5]]),cex = 1.5, col = "black") # Display Chart Title
    abline(h = par("usr")[3], col = "black") # Seal the black line at the bottom
    abline(v = par("usr")[2], col = "black") # Seal the black line on the right side
    
    #-----------------------------#
    # Canvas area 25: Draw legend #
    #-----------------------------#
    
    plot(0,0,col = "white",
         xlab = "",xaxt = "n", # Do not display the X-axis
         ylab = "",yaxt = "n") # Do not display the Y-axis
    
    legend("topleft",
           legend = c("Yes","No",
                      "Yes","No",
                      "Yes","No",
                      "Yes","No",
                      "Yes","No"),
           fill = c(c("#16a085","#e74c3c"),
                    c("#16a085","#e74c3c"),
                    c("#16a085","#e74c3c"),
                    c("#16a085","#e74c3c"),
                    c("#16a085","#e74c3c")),
           border = NA, # Legend color has no border
           bty = "n", # Legend has no border
           cex = 1.8,
           box.lwd = 3,
           x.intersp = 0.05,
           y.intersp = 1,
           text.width = 0.075, # Legend interval
           horiz = T) # Legend horizontal placement
    
    
  }
  table(ts$Risk)
  
}





