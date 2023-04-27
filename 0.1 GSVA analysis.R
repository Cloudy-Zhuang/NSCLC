# ## 3.PI3K/AKT/mTOR pathway GSVA analysis --------------------------------------------------
load("dat_input/NSCLC-exp-clin-data.Rdata")
mtor <- read.gmt("dat_input/HALLMARK_PI3K_AKT_MTOR_SIGNALING.v7.5.1.gmt")

NSCLC_fpkm <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)
#Poisson distribution was selected for sequencing data, and Gaussian distribution was selected for chip
mat <- as.matrix(NSCLC_fpkm)
mtor_list <- split(mtor$gene, mtor$term)

NSCLC.mtor.gsva.res <- gsva(expr=mat, 
                       mtor_list , 
                       kcdf="Poisson",
                       method = "gsva",
                       parallel.sz=10)


save(gsva_res,
     file="dat_output/NSCLC_mtor_gsva_res.Rdata")
