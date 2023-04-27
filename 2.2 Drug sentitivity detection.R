
library(pRRophetic)
library(tidyverse)
load("dat_input/药物预测准备数据.Rdata")
NSCLC_index_high %>% 
  select(which(substr(names(.),14,15)!="11"))->
  NSCLC_index_high.cancer

NSCLC_index_low %>% 
  select(which(substr(names(.),14,15)!="11"))->
  NSCLC_index_low.cancer


drug <- read.csv("dat_input/NSCLC drug list.csv") %>% 
  pull(x)


library(future)
library(future.apply)


data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
data(cgp2016ExprRma)
data(drugAndPhenoCgp)
possibleDrugs2016 <- unique( drugData2016$Drug.name)

drug_nsclc <- openxlsx::read.xlsx("dat_input/Pathway inhibitors information.xlsx",sheet = 2)


ind <- tolower(trimws(possibleDrugs2016))%in%tolower(trimws(drug$drug))

possibleDrugs2016[ind]
ind2 <- toupper(trimws(possibleDrugs2016))%in%toupper(trimws(drug$drug))
possibleDrugs2016[ind2]
plan(multisession, workers = 12)

predict_drug <- function(drug,dat){
  
  dat <- as.matrix(dat)
  
  predictedPtype=try(pRRopheticPredict(
    testMatrix=dat,
    drug=drug,
    tissueType = "all", 
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016"))
  
  if (class(predictedPtype)=="try-error") {
    return(NA)
  }else{return(predictedPtype)}
 
   gc()
  
}
drug <- c("Temsirolimus", "BEZ235",drug_nsclc)
drug_all_high<- lapply(drug,
                        predict_drug,
                        NSCLC_index_high.cancer)
drug_all_low <- lapply(drug,
                        predict_drug,
                        NSCLC_index_low.cancer)


