## WGCNA
load("dat_output/NSCLC_mtor_gsva_res.Rdata")
load("dat_output/NSCLC_gsva_res.Rdata")
load("dat_input/NSCLC-exp-clin-data.Rdata")
library(tidyverse)


# data
NSCLC_fpkm <- cbind(LUAD_fpkmanno,LUSC_fpkmanno)
NSCLC_clin <- rbind(LUAD_clin_tab,LUSC_clin_tab)

substr(rownames(NSCLC.mtor.gsva.res),14,15) ->ind1
substr(rownames(NSCLC.pi3k.gsva.res),14,15) ->ind2

  
NSCLC.mtor.gsva.res %>% 
  mutate(group=ifelse(ind1==11,"normal","cancer")) %>% 
  arrange(group)->
  pheno_mtor

NSCLC.pi3k.gsva.res %>% 
  mutate(group=ifelse(ind2==11,"normal","cancer")) %>% 
  arrange(group)->
  pheno_pi3k
# Merge gsva scoring matrix



# WGCNA process
library(WGCNA)

datExpr0 <- as.data.frame(t(NSCLC_fpkm))
test <- datExpr0[1:10,1:10]#Row is sample, column is gene
###################################################################################
### 1.Input conditions of goodSamplesGenes(Row is sample, column is gene)
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

### If it doesn't meet the standard, it needs to be screened

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    #printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


###Not enough, we can still do some screening work

#Turn around: row is gene, column is sample
raw_counts <- as.data.frame(t(datExpr0))

### 2.rm low counts data
low_count_mask <- rowSums(raw_counts) < ncol(raw_counts)
raw_counts <- raw_counts[!low_count_mask,]
sum(low_count_mask)


### 4.rm unchanged genes
# first, let's remove any genes with _zero_ variance since these are not
# going to help us, and may cause problems with some of the models
#Delete samples with expression variance ≤ 0
#Row is gene, column is sample
raw_counts <- raw_counts[apply(raw_counts, 1, var) > 0,]

### Cluster to see abnormal samples
sampleTree = hclust(dist(t(raw_counts)), method = "average")
### plot
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", 
     xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)



### Extract the intersection of expression quantity and trait information, and align the order of phenotypic data samples with the order of expression matrix
raw_counts %>% 
  select(rownames(pheno_mtor))->
  fpkm_dat_mtor
  
raw_counts %>% 
  select(rownames(pheno_pi3k))->
  fpkm_dat_pi3k

identical(fpkm_dat_mtor,fpkm_dat_pi3k)
#The two data are the same, just choose one



#The row is the sample and the column is the gene.
datExpr = as.data.frame(t(fpkm_dat_mtor))

test <- datExpr[1:10,1:10]
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(as.numeric(as.factor(pheno_mtor$HALLMARK_PI3K_AKT_MTOR_SIGNALING)),
                             signed = FALSE)
plotDendroAndColors(sampleTree2, 
                    traitColors,
                    groupLabels = names(pheno_mtor)[1], 
                    main = "Sample dendrogram and trait heatmap")

setwd("dat_output/")

save(datExpr,pheno_mtor,pheno_pi3k,file = "WGCNA_preparation_data.Rdata")

## WGCNA ANALYSIS
load("dat_output/WGCNA_preparation_data")
library(WGCNA)
### Set up multithreading
enableWGCNAThreads()
### Load the data of the previous step

test <- as.data.frame(datExpr)[1:20,1:20]
### Determination of soft threshold,
### Use the function pickSoftThreshold
### Give some values artificially
powers = c(c(1:10), seq(from = 12, to=30, by=2))
length(powers)
### The length of this power determines the number of rows of the final returned result
# Call the network topology analysis function
### Modules here are all used for expression data, and there is no trait data
### It takes a little time.
sft = pickSoftThreshold(datExpr, powerVector = powers,networkType="signed", verbose = 5)

test <- sft$fitIndices
### At this time, the result has been shown, but you can see it again by drawing
sft$powerEstimate
### Draw and render results

sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

### One-step network construction and module discovery # 
if (F) {
  ### power = 12, give the reason
  ##  maxBlockSize = 5000, It can be adjusted according to your computer, 16GB memory 20000, 32GB memory 30000
  ### TOMType 
  ### This step is the cornerstone of all analysis
  ### Prevent error reporting
  cor <- WGCNA::cor
  net = blockwiseModules(datExpr,
                         power = 12,
                         TOMType = "signed", 
                         minModuleSize = 30,
                         networkType="signed",
                         reassignThreshold = 0, 
                         mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "nsclc_wgcna", 
                         verbose = 3)
  cor<-stats::cor
  
  
  
  
  ### The returned net is a list, and the Module Epigene expression data is also in it
  table(net$colors)
  ### 0Represents genes that are not in the module
  
  ### draw
  # open a graphics window
  sizeGrWindow(12, 9)
  ### Label to color
  mergedColors = labels2colors(net$colors)
  table(mergedColors)
  # Plot the dendrogram and the module colors underneath
  plot(net$dendrograms[[1]])
  plotDendroAndColors(net$dendrograms[[1]], 
                      mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  ### Save data
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  geneTree = net$dendrograms[[1]]
  setwd("dat_output/")
  save(net,moduleLabels, moduleColors, geneTree, file = "NSCLC_module_detection.RData")
}


# Multistep method ---------------------------------------------------------------------

if (T) {
  
  softPower = 12 #NSCLCThe result is12
  adjacency = adjacency(datExpr, power = softPower)
  test <- adjacency[1:10,1:10]
  #  Code chunk 4
  #=====================================================================================
  
  # Turn adjacency into topological overlap
  TOM = TOMsimilarity(adjacency)
  class(TOM)
  dim(TOM)
  test <- TOM[1:10,1:10]
  dissTOM = 1-TOM
  
  #=====================================================================================
  #  Code chunk 5
  #=====================================================================================
  
  # Call the hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = "average")
  # Plot the resulting clustering tree (dendrogram)
  sizeGrWindow(12,9)
  plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04);
  
  #=====================================================================================
  #  Code chunk 6
  #=====================================================================================
  
  # We like large modules, so we set the minimum module size relatively high:
  minModuleSize = 30
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, 
                              distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  table(dynamicMods)
  
  #=====================================================================================
  #  Code chunk 7
  #=====================================================================================
  
  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  # Plot the dendrogram and colors underneath
  sizeGrWindow(8,6)
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  #=====================================================================================
  #  Code chunk 8
  #=====================================================================================
  ### Expression of similar module merging
  # Calculate eigengenes
  MEList = moduleEigengenes(datExpr, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs)
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average")
  # Plot the result
  sizeGrWindow(7, 6)
  plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
  
  #=====================================================================================
  #  Code chunk 9
  #=====================================================================================
  
  MEDissThres = 0.25
  # Plot the cut line into the dendrogram
  abline(h=MEDissThres, col = "red")
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr, 
                            dynamicColors, 
                            cutHeight = MEDissThres, 
                            verbose = 3)
  # The merged module colors
  mergedColors = merge$colors
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs
  
  #=====================================================================================
  #  Code chunk 10
  #=====================================================================================
  
  sizeGrWindow(12, 9)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  #=====================================================================================
  #  Code chunk 11
  #=====================================================================================
  
  # Rename to moduleColors
  moduleColors = mergedColors
  # Construct numerical labels corresponding to the colors
  ### why grey?
  colorOrder = c("grey", standardColors(50))
  moduleLabels = match(moduleColors, colorOrder)-1
  MEs = mergedMEs
  
  
  
setwd("dat_output/")
 
save(moduleLabels,
     mergedColors,
     geneTree,
     file = "module_step_by_step_detection.RData")
  
}





# Phenotype is associated with expression matrix --------------------------------------------------------------
load("dat_output/WGCNA_preparation_data.Rdata")
load("dat_output/NSCLC_mtor_gsva_res.Rdata")
load("dat_output/NSCLC_gsva_res.Rdata")
#load("dat_output/NSCLC_module_detection.RData")
load("dat_output/module_step_by_step_detection.RData")
library(tidyverse)
library(WGCNA)

substr(rownames(NSCLC.mtor.gsva.res),14,15) ->ind1
substr(rownames(NSCLC.pi3k.gsva.res),14,15) ->ind2

## Complete phenotypic data

NSCLC.mtor.gsva.res %>% 
  mutate(group=ifelse(ind1==11,"normal","cancer")) %>% 
  arrange(group)->
  pheno_mtor

NSCLC.pi3k.gsva.res %>% 
  mutate(group=ifelse(ind2==11,"normal","cancer")) %>% 
  arrange(group)->
  pheno_pi3k


# Recalculate MEs with color labels1
test <- moduleEigengenes(datExpr, moduleColors)

### Extract results
MEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
### Reordering, which changes the position of the column, has no practical significance, and is mainly used for drawing
MEs = orderMEs(MEs0)
#############################################################################
### Very important steps, very simple operation
### For relevance, what are the requirements of datTraits??
datTraits <- pheno_mtor
#datTraits <- pheno_pi3k
datTraits$condition <- factor(datTraits$group)

### Be careful... (I think we should pay attention to the order of samples)
design=model.matrix(~0+ datTraits$condition)
design <- as.data.frame(design)
colnames(design)=levels(datTraits$condition)
design
design$mtor <- as.numeric(datTraits$HALLMARK_PI3K_AKT_MTOR_SIGNALING)
#design$pi3k <- as.numeric(datTraits$WP_PI3KAKT_SIGNALING_PATHWAY)


### How to expand?
### GZ04_ TCGA tumor data mining - how to quantify tissue immune infiltration
### https://weidian.com/item.html?itemID=3636107475
#Module is related to phenotype
moduleTraitCor = cor(MEs, design, use = "p")
### calc p value
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
dd2 <- data.frame(table(mergedColors)) #
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
print(myp)
###The most important figure has been obtained
#######################################################################################
###Extract interested genes for mapping
###It must be extracted from the expression matrix
#module = "blue"
#Module="salmon" # salmon in pi3k akt path
#ModuleColors should be the color corresponding to each gene
module="turquoise"

moduleGenes = mergedColors==module
datExpr <- as.data.frame(datExpr)
dd <- datExpr[,moduleGenes]

test <- datExpr[1:20,1:20]
###############################################################################
### Explore the genes in this module

design
mylove = data.frame(design$mtor,design$cancer)
names(mylove) = c("mtor","cancer")

## select from 3rd
modNames = substring(names(MEs), 3)

### Calculate gene correlation: gene and module correlation
### module membership MM

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 
                                          nSamples))#NSamples is the sample size

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

### Gene Significance GS
###[mylove] The phenotype extracted from the phenotype matrix of my concern
###[geneTraitSignificance] is the expression matrix and interest phenotype correlation matrix
geneTraitSignificance = as.data.frame(cor(datExpr, mylove, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
                                          nSamples))

names(geneTraitSignificance) = paste("GS.", names(mylove), sep="")
names(GSPvalue) = paste("p.GS.", names(mylove), sep="");

moduleColors=mergedColors

### give example
#module = "blue"
#module = "salmon"
module="turquoise"
modNames = substring(names(MEs), 3)
column = match(module, modNames)
table(moduleColors)
unique(moduleColors)
### Extract the gene in the module
moduleGenes = moduleColors==module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
###It is essentially a scatter diagram
#Extract the correlation coefficient of interest module genes after all gene modules are related to the expression matrix
MM <- geneModuleMembership[moduleGenes, column]
#Extract the coefficient after correlation between expression matrix and interest phenotype
GS <- geneTraitSignificance[moduleGenes, 1]
x11()
verboseScatterplot(MM,GS,
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for PI3K/AKT Pathway",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

###############################################################################
### look for hub gene
## mtor
if (F) {
  # mtor
  MM  <- geneModuleMembership[moduleGenes, column]
  MMP <- MMPvalue[moduleGenes, column]
  GS <-  geneTraitSignificance[moduleGenes, 1]
  GSP <- GSPvalue[moduleGenes, 1]
  
  mydata_mtor <- data.frame(moduleGenes=colnames(datExpr)[moduleGenes],MM,MMP,GS,GSP)
  
  ## cancer
  MM  <- geneModuleMembership[moduleGenes, column]
  MMP <- MMPvalue[moduleGenes, column]
  GS <-  geneTraitSignificance[moduleGenes, 2]
  GSP <- GSPvalue[moduleGenes, 2]
  
  mydata_mtor_cancer <- data.frame(moduleGenes=colnames(datExpr)[moduleGenes],MM,MMP,GS,GSP)
  
  setwd("dat_output/")
  
  save( mydata_mtor,mydata_mtor_cancer,file = "Preliminary_results_of_blue_module_core_gene.Rdata")

mydata_mtor_cancer %>% 
  filter(abs(MM)> 0.8 & abs(GS)> 0.2)
  
  
  
}

# #Screening significant genes
# library(dplyr)
# mydata <- mydata %>% 
#   filter(MMP < 0.05,GSP < 0.05) %>% 
#   arrange(desc(MM),desc(GS)
rownames(NSCLC.pi3k.gsva.res) %>% 
  substr(14,15) ->sample_ind

NSCLC.pi3k.gsva.res %>% 
  mutate(group=ifelse(sample_ind=="11","normal","cancer")) %>% 
  ggplot(aes_string("group",names(.)[1]))+
  geom_boxplot()+
  ggpubr::stat_compare_means()

NSCLC.mtor.gsva.res %>% 
  mutate(group=ifelse(sample_ind=="11","normal","cancer")) %>% 
  ggplot(aes_string("group",names(.)[1],fill="group"))+
  geom_boxplot()+
  ggpubr::stat_compare_means(label.x = 1.2,
                             size=5)+
  ggsci::scale_fill_aaas()+
  theme_bw(base_size = 18)






