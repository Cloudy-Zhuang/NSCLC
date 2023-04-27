# Phenotype is associated with expression matrix --------------------------------------------------------------
load("dat_output/WGCNA_preparation_data.Rdata")
load("dat_output/NSCLC_mtor_gsva_res.Rdata")
load("dat_output/NSCLC_mtor_gsva_res.Rdata")
load("dat_output/module_step_by_step_detection.RData")




library(tidyverse)
library(WGCNA)

substr(rownames(NSCLC.mtor.gsva.res),14,15) ->ind1


## Complete phenotypic data

NSCLC.mtor.gsva.res %>% 
  mutate(group=ifelse(ind1==11,"normal","cancer")) %>% 
  arrange(group)->
  pheno_mtor
# Recalculate MEs with color labels1
test <- moduleEigengenes(datExpr, mergedColors)

### Extract results
MEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
### Reordering, which changes the position of the column, has no practical significance, and is mainly used for drawing
MEs = orderMEs(MEs0)
#############################################################################
### Very important steps, very simple operation
### For relevance, what are the requirements of datTraits?
# datTraits <- pheno_mtor
datTraits <- pheno_pi3k
datTraits$condition <- factor(datTraits$group)

### Be careful... 
design=model.matrix(~0+ datTraits$condition)
design <- as.data.frame(design)
colnames(design)=levels(datTraits$condition)
design

## This step is modified manually
design$mtor <- as.numeric(datTraits$HALLMARK_PI3K_AKT_MTOR_SIGNALING)

# Module is related to phenotype
moduleTraitCor = cor(MEs, design, use = "p")
### calc P-value
nSamples <- nrow(datExpr)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
### Magic operation!!
textMatrix =  paste(signif(moduleTraitCor, 2), "(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

### How to display the number of genes
dd1 <- data.frame(ME=names(MEs),color=substring(names(MEs),3))
dd2 <- data.frame(table(mergedColors))
colnames(dd2) <- c("color","num")
## merge
dd3 <- merge(dd1,dd2,by="color")
rownames(dd3) <- dd3$ME
dd3 <- dd3[names(MEs),]

### c(bottom, left, top, right)
par(mar = c(4, 12, 1, 1))
ynumbers <- paste0(names(MEs),paste0("(",dd3$num,")"))
ynumbers
myp <- labeledHeatmap(Matrix = moduleTraitCor,
                      xLabels = colnames(design),
                      yLabels = names(MEs),
                      ySymbols = ynumbers,
                      colorLabels = FALSE,
                      colors = blueWhiteRed(50),
                      textMatrix = textMatrix,
                      setStdMargins = FALSE,
                      cex.text = 0.8,
                      zlim = c(-1,1),
                      main = paste("Module-trait relationships"))
x11()

###The most important figure has been obtained
#######################################################################################
###Extract interested genes for mapping
###It must be extracted from the expression matrix

module="turquoise"
#module="tan"
moduleGenes = mergedColors==module
datExpr <- as.data.frame(datExpr)
dd <- datExpr[,moduleGenes]

test <- datExpr[1:20,1:20]
###############################################################################
### Explore the genes in this module

design
## Manual modification
mylove = data.frame(design$mtor,design$cancer)
names(mylove) = c("mtor","cancer")


## Select from the third place
modNames = substring(names(MEs), 3)

### Calculate gene correlation: gene and module correlation
### module membership MM

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 
                                          nSamples))#NSamples is the sample size

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

### Gene Significance GS
### The phenotypes extracted by mylove from the phenotypic matrix of its own concern
### GeneTraitSignificance is the expression matrix and interest phenotype correlation matrix
geneTraitSignificance = as.data.frame(cor(datExpr, mylove, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
                                          nSamples))

names(geneTraitSignificance) = paste("GS.", names(mylove), sep="")
names(GSPvalue) = paste("p.GS.", names(mylove), sep="");





module="turquoise"
modNames = substring(names(MEs), 3)
column = match(module, modNames)
table(mergedColors)
unique(mergedColors)
### Extract the gene in the module
moduleGenes = mergedColors==module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
###It is essentially a scatter diagram
#Extract the correlation coefficient of interest module genes after all gene modules are related to the expression matrix
MM <- geneModuleMembership[moduleGenes, column]
#Extract the coefficient after correlation between expression matrix and interest phenotype
GS <- geneTraitSignificance[moduleGenes, 2] #1 is pathway, 2 is tumor
x11()
mylabs='Gene significance for Cancer'
verboseScatterplot(MM,GS,
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = mylabs,
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

###############################################################################
### Find the hub gene by yourself
## mtor
if (F) {
  # mtor
  MM  <- geneModuleMembership[moduleGenes, column]
  MMP <- MMPvalue[moduleGenes, column]
  GS <-  geneTraitSignificance[moduleGenes, 1]
  GSP <- GSPvalue[moduleGenes, 1]
  
  merge_mtor <- data.frame(moduleGenes=colnames(datExpr)[moduleGenes],MM,MMP,GS,GSP)
  
  ## cancer
  MM  <- geneModuleMembership[moduleGenes, column]
  MMP <- MMPvalue[moduleGenes, column]
  GS <-  geneTraitSignificance[moduleGenes, 2]
  GSP <- GSPvalue[moduleGenes, 2]
  
  merge_mtor_cancer <- data.frame(moduleGenes=colnames(datExpr)[moduleGenes],MM,MMP,GS,GSP)
  
  setwd("dat_output")
  
  save( mydata_mtor,mydata_mtor_cancer,file = "Preliminary results of core genes in pi3k_akt_mtor_turquoise module after module merging.Rdata")
  
}


