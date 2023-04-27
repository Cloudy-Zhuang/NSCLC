#ROC diagram
#Roc calculation
rm(list = ls())
library(tidyverse)
library(pROC)
library(plotROC)
library(wesanderson)
load("dat_input/Data prepared for roc.Rdata")


# Create function --------------------------------------------------------------------
if (T) {
  roc_count <- function(gene,outcome,data){
    
    
    roc.data = roc(data[,outcome],data[,gene],
                   percent=T,plot=F, grid=T,lty=1,
                   print.auc=F)
    #P-value extraction
    se <- sqrt(var(roc.data))  # Obtain SE of AUC
    b <- roc.data$auc - .5
    z <- (b / se)  # Calculate Z value
    p.value=2 * pt(-abs(z), df=Inf)  ## two-sided tes
    auc <- round(as.numeric(roc.data$auc)/100,2)
    #Other parameter extraction
    index <-coords(roc.data,"best")
    Threshold=index[1,1]
    Specificity=index[1,2]/100
    Sensitivity=index[1,3]/100
    Youden_index = Sensitivity + Specificity -1
    data.frame(gene,auc, Specificity,Sensitivity,Youden_index,p.value)
  }
  
  
  
  
  
  
  
  roc_plot <- function(gene,outcome,data,write.out=F,color="red"){
    
    if (write.out!=F) {
      pdf(paste0(gene,".pdf"),height = 4,width = 4)
    }
    
    par(mar=c ( 5, 4, 4, 2)-1)
    roc.data = roc(data[,outcome],data[,gene],
                   percent=T,plot=T, grid=T,lty=1,
                   print.auc=F,col=color,
                   thresholds="best",
                   print.thres="best")
    #P-value extraction
    se <- sqrt(var(roc.data))  # Obtain SE of AUC
    b <- roc.data$auc - .5
    z <- (b / se)  # Calculate Z value
    p.value=2 * pt(-abs(z), df=Inf)  ## two-sided tes
    p=ifelse(p.value>0.05,"p>0.05",
             ifelse(p.value<0.001,"p<0.001",
                    ifelse(p<0.01,"p<0.01","p<0.05")))
    #Other parameter extraction
    index <-coords(roc.data,"best")
    Threshold=paste0("Threshold=",round(index[1],1))
    Specificity=paste0("Specificity=",round(index[2],1),"%")
    Sensitivity=paste0("Sensitivity=",round(index[3],1),"%")
    title(main=gene,cex=1,outer=F)
    text(25,25,
         paste0("AUC",":",round(as.numeric(roc.data$auc)/100,2),"\n",
                p,"\n",
                Threshold,"\n",
                Specificity,"\n",
                Sensitivity),
         col="black")
    
    if (write.out!=F) {
      dev.off() 
    }
    
    

    
  }
}

roc.pre.dat %>% 
  select(group,everything())->
  roc_data


roc.res <- map_df(names(roc_data[-1]),
       roc_count,
       data=roc_data,
       outcome = "group")

roc.res %>% 
  mutate(across(where(is.numeric),round,3)) %>% 
  DT::datatable()


openxlsx::write.xlsx(roc.res,"dat_output/Important gene roc results.xlsx")

roc.res %>% 
  DT::datatable()

# draw figure----------------------------------------------------------------------


lapply(names(roc_data[-1]),
       roc_plot,
       data=roc_data,
       outcome = "group",
       write.out=T,
       color="#fb6501")

