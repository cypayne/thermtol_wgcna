###  Differential gene expression analysis with DESeq2
#     Need:
#           kallisto output (abundance.h5 file per sample)
#           GFF/GTF with annotated transcripts
#           sample file

## Tutorials: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/

library(tximportData)
library(tximport)
library(GenomicFeatures)
library(readr)
library(DESeq2)
library(rhdf5)

# collect kallisto output
setwd("~/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/kallisto_output/")

######  Process input data - DGE analysis with DESeq2 ######
# Want to pass in variance stabilized data from DESeq2 to WGCNA

### Pre-process kallisto read abundance files ### 

# XXX specify tissue
tissue <- "brain"

# XXX specify kallisto directory
dir <- "./kall2birch-posttrim-w-mito"

## Read in sample files
samples <- read.table("./TT_samples.txt", header = TRUE)

# make sure that non-continuous variables are cast as factors
samples$temp <- factor(samples$temp)
samples$species<-as.factor(samples$species)

# subset data by tissue of interest
samples <- samples[samples$tissue == tissue, ]
samples

# collect files from specified kallisto directory
files <- file.path(dir, paste(samples$file_basename,"kallisto",sep="_"), "abundance.h5")
names(files) <- paste0(samples$sample)
files

# make transcript annotation object with GFF/GTF
# note: this file should only include "transcript" elements, remove "CDS", "start", "stop", etc elements
txdb <- makeTxDbFromGFF(file="../refs/Xiphophorus_maculatus_LG.Xipmac4.4.2.81_unique-tx_w-mito.gtf", format="gtf")

# create tx2gene table, matching transcript name to geneid
k <- keys(txdb, keytype = "TXNAME")
tx2gene <-  AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# get gene-level counts, since kallisto only provides transcript-level data
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
names(txi)
head(txi$counts)

# ### Choosing whether to use an interaction factor in the design ###
# design <- as.formula(~ species + temp + species:temp)
# dds <- DESeqDataSetFromTximport(txi, colData = samples, design = design)
# dds <- DESeq(dds)
# redDesign <- as.formula(~ species + temp )
# # Compare the designs
# ddsComp <- DESeq(dds, test="LRT", reduced=redDesign)
# res.33CbirVmal <- lfcShrink(ddsComp, coef="species_bir_vs_mal", type="ashr")
# res.mal33V22 <- lfcShrink(ddsComp, coef="temp_33C_vs_22C", type="ashr")
# 
# resSig <- subset(res.33CbirVmal, padj < 0.05)
# resSig <- subset(res.mal33V22, padj < 0.05)
# # look at how many genes have sig difference
# dim(resSig)
# 
# ### done choosing design ###


### DGE with DESeq2 to get transformed counts ### 

# set design 
design <- as.formula(~ replicate + species + temp)

# instantiate the DESeqDataSet
dds <- DESeqDataSetFromTximport(txi,
                                colData = samples,
                                design = design)

# set the reference/control group, to be compared against
dds$species <- relevel(dds$species, ref = "Xmal")
dds$temp <- relevel(dds$temp, ref = "22.5")

# run DESeq for size factor est, dispersion est, and 
# negative binomial GLM fitting and significance test
dds <- DESeq(dds)

# look at names of estimated effects
resultsNames(dds)

#variance stabilized data for WGCNA
#vsd = getVarianceStabilizedData(dds)
vsd <- assay(vst(dds, blind=FALSE))

write.csv(vsd,"BrainTT-species.temp_vsd.csv",quote=F)

#add log 2 fold
logvsd = log2(vsd)

# output logvsd
write.csv(logvsd,"BrainTT-species.temp_log2vsd.csv",quote=F)


###### Running WGCNA ######

# Load the package
library(WGCNA)
options(stringsAsFactors = FALSE)

# load vsd counts
vsd <- read.csv("BrainTT-species.temp_vsd.csv",header=T)
rownames(vsd) <- vsd$X
vsd <- vsd[,-1]

#quick look at data
dim(vsd)

# transpose, put into dataframe
datExpr0 = as.data.frame(t(vsd))


### Filter genes and samples ###

# check for missing entries, weights below a threshold, 
# and zero-variance genes, returns list of genes/samples
# that pass these filters
gsg = goodSamplesGenes(datExpr0, minNSamples = 17, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average");

## plot sample tree to look for outliers
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

## Note: 
#S9_TSR1-22C-malxbirchF1-brain_9-C4 looks like it could be an outlier

## Remove sample outliers by choosing height for branch cut

# Plot a line to show the cutoff
abline(h = 15, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==0) # change the cluster value, in this case only 1 cluster=0
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


### Make sample dendogram with Trait heat map ###

## Prepare traits

# XXX To one-hot encode your traits:
library(data.table)
Traits <-samples[, c("sample","species","temp")]
one_hotTraits <- one_hot(as.data.table(Traits), cols =c("species", "temp"))
allTraits<-one_hotTraits

# XXX Otherwise:
allTraits<-samples[, c("sample","species","temp")]
#rownames(allTraits)<-samples$sample

expSamples = rownames(datExpr);
traitRows = match(expSamples, allTraits$sample);
#traitRows = data.frame(expSamples, allTraits$sample)

datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];

## Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
#traitColors = numbers2colors(cbind(as.factor(datTraits[,1]),as.factor(datTraits[,2])), signed = FALSE);
traitColors = numbers2colors(cbind(as.factor(datTraits[,1]),as.factor(datTraits[,2]),as.factor(datTraits[,3]),as.factor(datTraits[,4]),as.factor(datTraits[,5])), signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
pdf(file='TT-brain-WGCNA_sample-tree.pdf',height=7,width=9.5)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

## Save objects for subsequent easy loading:
save(datExpr, datTraits, file = "BrainThermalData-01-dataInput.RData") 

lnames = load(file = "BrainThermalData-01-dataInput.RData");


### Choose Soft-thresholding power ###

## Test a set of soft-thresholding powers
## choose lowest power for which scale-free topoplogy index (SFT.R.sq) reaches 0.9 or inflection point of plot
## https://support.bioconductor.org/p/87024/
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

## Plot the results:
pdf(file='TT-brain-WGCNA_pickSoftThreshold.pdf',height=7,width=9.5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Plot 1: Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Plot 2: Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

## Looks like 7 is the asymptoting soft-threshold

### Cluster genes into modules ###

## We'll use a single block (i.e. all 19k genes in one go, use maxBlockSize=20000)
# maxBlockSize = 20000, blockSizePenaltyPower = Inf : if maxBlockSize > total tx
# maxBlockSize = 10000, blockSizePenaltyPower = 5 : otherwise
net = blockwiseModules(datExpr, power = 7,
                       TOMType = "unsigned", minModuleSize = 20,
                       maxBlockSize = 20000, blockSizePenaltyPower = Inf,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = FALSE, pamRespectsDendro = FALSE,
                       loadTOMs = FALSE,
                       saveTOMFileBase = "TT-Brain-TOM",
                       verbose = 3)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath
pdf(file='TT-brain-WGCNA_blockwiseModule-dendrogram_onehot.pdf',height=8.5,width=12)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#save modules as an Rdata object
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "BrainThermal-02-networkConstruction-auto.RData")


### Correlate modules with traits of interest ###

# Load the expression and trait data saved in the first part
lnames = load(file = "BrainThermalData-01-dataInput.RData");
lnames
# Load network data saved in the second part
lnames = load(file = "BrainThermal-02-networkConstruction-auto.RData");
lnames

## Get module eigenvalues and calculate correlation with each trait
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# if not using one-hot encoded traits
moduleTraitCor = cor(MEs, cbind(datTraits$species,datTraits$temp), use = "p");
# if using one-hot encoded traits
moduleTraitCor = cor(MEs, cbind(datTraits$species_malxbirchF1,datTraits$species_Xbirch,datTraits$species_Xmal,datTraits$temp_22.5,datTraits$temp_33.5), use = "p");
#moduleTraitCor = cor(MEs, cbind(datTraits$temp_22.5,datTraits$temp_33.5), use = "p");

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

## Display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf(file='TT-brain-WGCNA_module-trait-heatmap.pdf',height=10,width=12)
par(mar = c(6, 8.5, 3, 3))
# Plot heatmap of ME-trait correlations
hmcolors <- colorRampPalette(c("deepskyblue3", "white", "gold"))(n = 100)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = hmcolors,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


### Subset MEs with significant trait correlations ###
## For one-hot encoded traits

# Subset species related modules:
subset(moduleTraitPvalue,moduleTraitPvalue[,1]<0.05) # f1
subset(moduleTraitPvalue,moduleTraitPvalue[,2]<0.05) # xbirch
subset(moduleTraitPvalue,moduleTraitPvalue[,3]<0.05) # xmal
# Subset temp related modules:
subset(moduleTraitPvalue,moduleTraitPvalue[,4]<0.05) #22.5c
subset(moduleTraitPvalue,moduleTraitPvalue[,5]<0.05) #33.5c

## Output significant ME pvalues by trait
colnames(moduleTraitPvalue) <- c("F1","xbirch","xmal","temp22.5c","temp33.5")
write.csv(moduleTraitPvalue,"TT-brain-WGCNA_MEtraitpvals.csv")
write.csv(subset(moduleTraitPvalue,moduleTraitPvalue[,1]<0.05),"TT-brain-WGCNA_MEtraitpvals_sig-f1.csv",quote=F)
write.csv(subset(moduleTraitPvalue,moduleTraitPvalue[,2]<0.05),"TT-brain-WGCNA_MEtraitpvals_sig-xbirch.csv",quote=F)
write.csv(subset(moduleTraitPvalue,moduleTraitPvalue[,3]<0.05),"TT-brain-WGCNA_MEtraitpvals_sig-xmal.csv",quote=F)
write.csv(subset(moduleTraitPvalue,moduleTraitPvalue[,4]<0.05),"TT-brain-WGCNA_MEtraitpvals_sig-temp.csv",quote=F)

## Plot heatmap for modules with at least one significant trait relationship
sigmoduleTraitCor <- subset(moduleTraitCor,moduleTraitPvalue[,1] <0.05 | moduleTraitPvalue[,2] <0.05 | moduleTraitPvalue[,3] <0.05 | moduleTraitPvalue[,4] <0.05)
sigmoduleTraitPvalue<- subset(moduleTraitPvalue,moduleTraitPvalue[,1] <0.05 | moduleTraitPvalue[,2] <0.05 | moduleTraitPvalue[,3] <0.05 | moduleTraitPvalue[,4] <0.05)
sigMEs <-subset(names(MEs),moduleTraitPvalue[,1] <0.05 | moduleTraitPvalue[,2] <0.05 | moduleTraitPvalue[,3] <0.05 | moduleTraitPvalue[,4] <0.05)

textMatrix =  paste(signif(sigmoduleTraitCor, 2), sep = "");
dim(textMatrix) = dim(sigmoduleTraitCor)
hmcolors <- colorRampPalette(c("deepskyblue3", "white", "gold"))(n = 100)
# Plot heatmap
pdf(file='TT-brain-WGCNA_module-trait-heatmap_sigMEs.pdf',height=7.5,width=5.5)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = sigmoduleTraitCor,
               xLabels = names(datTraits),
               yLabels = sigMEs,
               ySymbols = sigMEs,
               colorLabels = FALSE,
#               colors = blueWhiteRed(50),
               colors = hmcolors,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


###### Output module membership for all modules with significant trait relationship ######

### Get summary output for each significant module ###

## make outfiles for all significant modules, in bulk:
header="TT-brain-WGCNA-onehot_"
for( MEname in sigMEs ) { 
  MEname <- substring(MEname, 3)
  print(MEname)
  ME<-names(datExpr)[moduleColors==MEname]
  write.csv(ME,paste0(header,MEname,"_genes.csv"))
}

# save sigMEs object for future use
save(sigMEs, file = "TT-brain-WGCNA_sigMEs.RData") 

## make outfiles for one module at a time:
temp_MEsalmon<-names(datExpr)[moduleColors=="salmon"]
write.csv(temp_MEsalmon,"BrainTT-temp_MEsalmon-0.78_genes.csv")

### Plot correlation between module membership and trait ###

## get module membership and p values for each gene
# simple Pearson's correlation
# yields same results as datKME <- signedKME(datExpr, MEs, outputColumnName="MM.")
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
write.csv(geneModuleMembership,paste0(header,"geneModuleMembership-temp33.5.csv"))
write.csv(MMPvalue,paste0(header,"geneModuleMembership_pvals-temp33.5.csv"))

## For temperature:
temp = as.data.frame(datTraits$temp_33.5);
names(temp) = "temp33.5"
modNames = substring(names(MEs), 3)

# get gene-temperature correlation and p values for each gene
geneTraitSignificance = as.data.frame(cor(datExpr, temp, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(temp), sep="");
names(GSPvalue) = paste("p.GS.", names(temp), sep="");
write.csv(geneTraitSignificance,paste0(header,"geneTraitSig-temp33.5.csv"))
write.csv(GSPvalue,paste0(header,"geneTraitSig_pvals-temp33.5.csv"))

# Plot module membership vs gene signficance for temp trait for individual modules
module = "palevioletred3"
column = match(module, modNames)
moduleGenes = moduleColors==module
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for temp",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# Plot kWithin (intramodular) connectivity and gene-trait significance for significant modules
colorlevels=unique(sigMEs)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorlevels))){
  whichmodule=colorlevels[[i]];
  restrict1 = (sigMEs==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     geneTraitSignificance$temp33.5[restrict1], col=sigMEs[restrict1],
                     main=whichmodule,
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}



## For species: 
species = as.data.frame(cbind(datTraits$species_malxbirchF1,datTraits$species_Xbirch,datTraits$species_Xmal))
names(species) = c("f1","xbirch","xmal")
modNames = substring(names(MEs), 3)

# get gene-species correlation and p values for each gene
geneTraitSignificance = as.data.frame(cor(datExpr, species, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(species), sep="")
names(GSPvalue) = paste("p.GS.", names(species), sep="")
write.csv(geneTraitSignificance,paste0(header,"geneTraitSig-onehot-species.csv"))
write.csv(GSPvalue,paste0(header,"geneTraitSig_pvals-onehot-species.csv"))

# Plot module membership vs gene signficance for species trait for individual modules
module = "blue"
column = match(module, modNames)
moduleGenes = moduleColors==module

pdf(file='TT-brain-WGCNA_0.95-palevioletred3_cor-plot.pdf',height=5,width=5)
par(mfrow = c(1,1))
# XXX change the geneTraitSignificance column to be plotted
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for species",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()


### Look at connectivity and significant genes in these modules ###

## Identify single top hub gene in a module
chooseTopHubInEachModule(datExpr,colorh="MEsalmon") # ENSXMAG00000011516

## Identify genes with high significance and high intramodular connectivity in a module (i.e. hub genes)
hub.genes = abs(geneTraitSignificance$GS.f1)> 0.2 & abs(geneModuleMembership$MEsalmon)> 0.85
table(hub.genes)
dimnames(data.frame(datExpr))[[2]][hub.genes]

## Get connectivity scores (from: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-07-Membership.pdf)
# the following is the same as running moduleCon<-intramodularConnectivity.fromExpr(datExpr,colors=moduleColors,corFnc = "cor", corOptions = "use='p'",networkType = "unsigned", power = 6)
connect=abs(cor(datExpr,use="p"))^6
connectivity_data=intramodularConnectivity(connect, moduleColors)
head(connectivity_data) # Higher values = more connectivity
dim(connectivity_data)
write.csv(connectivity_data,paste(header,"intramod-connectivity_genes.csv"))

## Get hub genes
# Two methods were used by Morgan et al 2020 to define hub genes:
#   1. any gene with module membership >= 0.85 (correlation between expression of gene and the module eigenvalue, Horvath & Dong 2008)
#   2. per module, top 5 most connected genes based on 
#       average number of neighbors/avg connectivity of nodes and network density/overall module connectivity
#       authors used Cytoscape to do this, as in Shannon et al 2003
#       use intramodular connectivity (kIN -- high module membership usually means high kIN

## kTotal (connectivity with all genes) and kWithin (connectivity with genes within module) quantiles have same number (as expected) but with some overlapping, some different members
kTotal_q0.99 <- subset(connectivity_data,connectivity_data$kTotal>quantile(connectivity_data$kTotal,0.99)) #190
kWithin_q0.99 <-subset(connectivity_data,connectivity_data$kWithin>quantile(connectivity_data$kWithin,0.99)) #190
kTotal_q0.95 <- subset(connectivity_data,connectivity_data$kTotal>quantile(connectivity_data$kTotal,0.95)) #946
kWithin_q0.95 <-subset(connectivity_data,connectivity_data$kWithin>quantile(connectivity_data$kWithin,0.95)) #946

#highconnectivity<-subset(row.names(connectivity_data),connectivity_data$kDiff>quantile(connectivity_data$kDiff,0.99))
#subset(Alldegrees1,rownames(Alldegrees1)=="ENSXMAG00000000428")

## Subset all genes with module membership values of +/-0.9 or greater
MM0.9 <- geneModuleMembership[rowSums(abs(geneModuleMembership) >= 0.9) >= 1, ] # 1036
MM0.85 <- geneModuleMembership[rowSums(abs(geneModuleMembership) >= 0.85) >= 1, ] # 2778

## Designate hub genes as genes in both kTotal_q0.95 and MM0.9 subsets
hub.genes <- intersect(row.names(kTotal_q0.95),row.names(MM0.9)) # 311
hub.genes <- intersect(row.names(kWithin_q0.95),row.names(MM0.9)) # 228
hub.genes <- intersect(row.names(kTotal_q0.95),row.names(MM0.85)) # 705
write.csv(connectivity_data[hub.genes,],paste0(header,"hubgenes_kTotal0.95_MM0.85.csv"))
hub.genes <- intersect(row.names(kWithin_q0.95),row.names(MM0.85)) # 600
write.csv(connectivity_data[hub.genes,],paste0(header,"hubgenes_kWithin0.95_MM0.85.csv"))

#hub.genes <- row.names(kWithin_q0.95)

## Identify which hub genes are found in modules of interest
MEsalmon_genes<-names(datExpr)[moduleColors=="salmon"]
intersect(hub.genes,MEsalmon_genes)
MEgrey60_genes<-names(datExpr)[moduleColors=="grey60"]
intersect(hub.genes,MEgrey60_genes)
greenyellowConn<-subset(hub.genes,X %in% MEgreenyellow_genes$x)

