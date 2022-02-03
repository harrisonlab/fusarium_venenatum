# Coexpression analysis

This analysis was done without collapsing replicas

I did this on my local machine in the WGCNA_VP_analysis

```r
library(WGCNA)
library(wesanderson) 

outdir <-setwd("./")

options(stringsAsFactors = FALSE)
# Input data
data=read.csv("iUK_data_vst_F.txt", sep="\t")
dim(data)
names(data)
# Transpose
datExpr0 = as.data.frame(t(data[, c(2:121)]))
names(datExpr0) = data$ID
rownames(datExpr0) = names(data)[c(2:121)]
# Check missing values
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

# Remove any genes and samples that do not pass the cut

if (!gsg$allOK){
    # Print items removed to file
    if (sum(!gsg$goodGenes) > 0){
        genes_to_remove <- (paste("Removing genes:", paste(names(datExpr0)[
        !gsg$goodGenes], collapse = "\n")))
        gfile <- paste(outdir, "removed_genes.txt", sep = "/")
        write(genes_to_remove, file = gfile)
    }
    if (sum(!gsg$goodSamples) > 0){
        samples_to_remove <- (paste("Removing samples:",
        paste(rownames(datexpr0)[!gsg$goodSamples], collapse = "\n")))
        sfile <- paste(outdir, "removed_samples.txt", sep = "/")
        write(samples_to_remove, file = sfile)
    }
    # Remove items that fail QC
    datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Cluster to check outliers

sampletree <- hclust(dist(datExpr0), method = "average")
file <- paste(outdir, "sample_clustering.pdf", sep = "/")
pdf(file, width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampletree, main = "Sample clustering to detect outliers", sub = "",
xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# Remove outlier samples, the height may need changing so be sure to check
abline(h = 30000, col = "red")
dev.off()

# Determine clusters under the line.
clust <- cutreeStatic(sampletree, cutHeight = 300000, minSize = 10)
table(clust)

# Clust 1 contains the samples we want to keep
keepsamples_clust1 <- (clust == 1)
datExpr1 <- datExpr0[keepsamples_clust1, ]
ngenes <- ncol(datExpr1)
nsamples <- nrow(datExpr1)

# Save a table with the expression data of all the samples
Rfile <- paste(outdir, "Cleaned_data.RData", sep = "/")
save(datExpr0, file = Rfile)

# Testing soft-thresholding power values

powers <- c(c(1:40))
#sft <- pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
sft <- pickSoftThreshold(datExpr0, powerVector = powers, networkType = "signed", verbose = 5,RsquaredCut=0.96)

# Draw a plot to allow manual picking of sft value

cex1 <- 0.9
file <- paste(outdir, "sft_testing.pdf", sep = "/")
pdf(file, height = 9, width = 12)
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
xlab = "Soft Threshold (power)",
ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)",
ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1,
col = "red")
dev.off()

write.csv(sft, file = "sft_indices.csv")

# here we define the adjacency matrix using soft thresholding with beta=6
ADJ1=abs(cor(datExpr0,use="p"))^14
k=softConnectivity(datExpr=datExpr0,power=14,type = "signed", verbose = 5)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(5,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check Scale-free Topology\n")

#adjacency <- adjacency(datExpr0, power = 6)
adjacency <- adjacency(datExpr0, type = "signed", power = 14)

# Topological Overlap Matrix (TOM)

# I test a few TOMtype parameters
tom <- TOMsimilarity(adjacency,TOMType = "signed")
# tom2 <- TOMsimilarity(adjacency) 
# tom3 <- TOMsimilarity(adjacency,TOMType = "unsigned")
# tom4 <- TOMsimilarity(adjacency,TOMType = "signed Nowick")

# The previous 4 gave the same results. This one is different but it is experimental and should be ignored
# tom5 <- TOMsimilarity(adjacency,TOMType = "signed 2")

disstom <- 1 - tom

# Clustering using TOM

#dist() will actually calculate the distance matrix whereas as.dist() only will try to coerce a object to a distance matrix.

#genetree2 <- hclust(as.dist(disstom), method = "average")
genetree <- hclust(dist(disstom), method = "average")
genetree3 <- hclust(d=as.dist(disstom), method = "average")

file <- paste(outdir, "clustering_tree.pdf", sep = "/")
pdf(file, height = 9, width = 12)
sizeGrWindow(12, 9)
plot(genetree2, xlab = "", sub = "", main = "Gene clustering on TOM-based
dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

# Cut clustering tree into several modules

# This gives module 0 for unassigned.
dynamicmods <- cutreeDynamic(dendro = genetree2, distM = disstom, deepSplit = 2, pamStage = TRUE,
pamRespectsDendro = TRUE, minClusterSize = 30)
table(dynamicmods)

# dynamicmods <- cutreeDynamic(dendro = genetree2, distM = disstom, deepSplit = 2,
# pamRespectsDendro = FALSE, minClusterSize = 30)
# table(dynamicmods)

power = 14
net = blockwiseModules(datExpr0, power = 14, checkMissingData = FALSE,
TOMType = "signed", minModuleSize = 30,
reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,
saveTOMFileBase = "Lei_TOM",
verbose = 3)
table(net$colors) 


# Plot modules on clustering tree, allows sanity check of min_mod_size value

file <- paste(outdir, "clustering_tree_with_modules.pdf", sep = "/")
pdf(file, height = 9, width = 12)
dynamiccolours <- labels2colors(dynamicmods)
table(dynamiccolours)
plotDendroAndColors(genetree2, dynamiccolours, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colours")
dev.off()

# Merging modules with similar expression patterns

melist <- moduleEigengenes(datExpr0, colors = dynamiccolours)
# No need to set softPower
# melist2 <- moduleEigengenes(datExpr0, colors = dynamiccolours,  softPower = 12)

mes <- melist$eigengenes
mediss <- 1 - cor(mes)
metree <- hclust(as.dist(mediss), method = "average")
file <- paste(outdir, "clustering_tree_with_merged_modules.pdf", sep = "/")
pdf(file, height = 9, width = 12)
plot(metree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h = 0.25, col = "red")
dev.off()

# Plot a comparison of merged and unmerged modules
merge <- mergeCloseModules(datExpr0, dynamiccolours, unassdColor="grey", cutHeight = 0.25,verbose = 3)
mergedcolours <- merge$colors
mergedmes <- merge$newMEs
file <- paste(outdir, "clustering_tree_compare_modules2.pdf", sep = "/")
pdf(file, height = 9, width = 12)
plotDendroAndColors(genetree2, cbind(dynamiccolours, mergedcolours),
c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

# Save output for further analyses

modulecolours <- mergedcolours
colourorder <- c("grey", standardColors(50))
modulelabels <- match(modulecolours, colourorder) - 1
mes <- mergedmes
Rfile <- paste(outdir, "modules.RData", sep = "/")
save(mes, modulelabels, modulecolours, dynamiccolours, genetree, tom, file = Rfile)


# Load list of transcript IDs and write out merged modules

transcripts <- names(datExpr0)
for (module in unique(modulecolours)){
    modgenes <- (modulecolours == module)
    modtranscripts <- transcripts[modgenes]
    filename <- paste("Genes_in_", module, ".txt", sep = "")
    mergedir <- paste(outdir, "merged_modules", sep = "/")
    file <- paste(mergedir, filename, sep = "/")
    dir.create(mergedir)
    write.table(as.data.frame(modtranscripts), file = file, row.names = FALSE,
    col.names = FALSE)
}

# Write out unmerged modules

for (module in unique(dynamiccolours)){
        modgenes <- (dynamiccolours == module)
        modtranscripts <- transcripts[modgenes]
        filename <- paste("Genes_in_", module, ".txt", sep = "")
        unmergedir <- paste(outdir, "unmerged_modules", sep = "/")
        file <- paste(unmergedir, filename, sep = "/")
        dir.create(unmergedir)
        write.table(as.data.frame(modtranscripts), file = file,
        row.names = FALSE, col.names = FALSE)
    }


# Associate modules to treatments


library(WGCNA)

options(stringsAsFactors = FALSE)

# Load the expression and trait data saved in the first part
lnames = load(file = "Cleaned_data.RData")
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "modules.RData")
lnames

# Create a matrix with 1 or 0 defining Media type of each sample
colData<-read.csv("colData.txt", sep="\t")
dim(colData)
names(colData)

# remove columns that hold information we do not need.
allTraits = colData[, -c(2,3,4)]
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.
expData = rownames(datExpr0)
traitRows = match(expData, allTraits$Sample.name)
datTraits = allTraits[traitRows, -1]
rownames(datTraits) = allTraits[traitRows, 1]
collectGarbage()

# Define numbers of genes and samples
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

# Merged modules

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, modulecolours)$eigengenes
MEs = orderMEs(MEs0)

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(12,9)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# Edit to resize
par(mar = c(9, 16, 3, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = wes_palette("Zissou1", 20, type = "continuous"),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-Condition relationships"))

# For unmerged modules

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, dynamiccolours)$eigengenes
MEs = orderMEs(MEs0)

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(12,9)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# Edit to resize
par(mar = c(9, 16, 3, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = wes_palette("Zissou1", 20, type = "continuous"),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-Condition relationships"))
```

```r
# Plot expression
darkorange<-read.table("Genes_in_turquoise.txt",header=FALSE)[,1]
vstRep2<-varianceStabilizingTransformation(dds,blind=FALSE,fitType="parametric")
mat <- assay(vstRep2)[darkorange, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vstRep2)[c("Media")])
pal <- wes_palette("Zissou1", 10, type = "continuous")
Z<-pheatmap(mat,color=pal,annotation_col = anno)
```

# MYRO

```r
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$MYRO)
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr0, weight, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

#Summary output of network analysis results

names(datExpr0)

#names(datExpr0)[modulecolours=="green"]

annot = read.csv(file = "WT_minion_annotations_geneID_short.csv")
dim(annot)
names(annot)
probes = names(datExpr0)
probes2annot = match(probes, annot$transcript_id)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))

# Create the starting data frame with all the annotations, GS.weight and p.GS.weight
geneInfo0 = data.frame(transcript_id = probes,
swissprot_org = annot$swissprot_org[probes2annot],
swissprot_gene = annot$swissprot_gene[probes2annot],
swissprot_function = annot$swissprot_function[probes2annot],
interproscan = annot$interproscan[probes2annot],
TFs = annot$TFs[probes2annot],
Antismash_type = annot$Antismash_type[probes2annot],
Antismash_cluster = annot$Antismash_cluster[probes2annot],
moduleColor = modulecolours,
geneTraitSignificance,
GSPvalue)

# Order modules by their significance for weight (from smallest to largets p-val). I uses MYRO weight
modOrder = order(-abs(cor(MEs, weight, use = "p")))
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
oldNames = names(geneInfo0)
geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
MMPvalue[, modOrder[mod]]);
names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight))
geneInfo = geneInfo0[geneOrder, ]
#This data frame can be written into a text-format spreadsheet, for example by
write.csv(geneInfo, file = "geneInfo_MYRO.csv")
```

# MYRO_KO

```r
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$MYRO_tri5_KO)
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr0, weight, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

#Summary output of network analysis results

names(datExpr0)

#names(datExpr0)[modulecolours=="green"]

annot = read.csv(file = "WT_minion_annotations_geneID_short.csv")
dim(annot)
names(annot)
probes = names(datExpr0)
probes2annot = match(probes, annot$transcript_id)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))

# Create the starting data frame with all the annotations, GS.weight and p.GS.weight
geneInfo0 = data.frame(transcript_id = probes,
swissprot_org = annot$swissprot_org[probes2annot],
swissprot_gene = annot$swissprot_gene[probes2annot],
swissprot_function = annot$swissprot_function[probes2annot],
interproscan = annot$interproscan[probes2annot],
TFs = annot$TFs[probes2annot],
Antismash_type = annot$Antismash_type[probes2annot],
Antismash_cluster = annot$Antismash_cluster[probes2annot],
moduleColor = modulecolours,
geneTraitSignificance,
GSPvalue)

# Order modules by their significance for weight (from smallest to largets p-val). I uses MYRO weight
modOrder = order(-abs(cor(MEs, weight, use = "p")))
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
oldNames = names(geneInfo0)
geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
MMPvalue[, modOrder[mod]]);
names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight))
geneInfo = geneInfo0[geneOrder, ]
#This data frame can be written into a text-format spreadsheet, for example by
write.csv(geneInfo, file = "geneInfo_MYRO_KO.csv")
```

# Beet_derived_thick_juice

```r
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$Beet_derived_thick_juice)
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr0, weight, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

#Summary output of network analysis results

names(datExpr0)

#names(datExpr0)[modulecolours=="green"]

annot = read.csv(file = "WT_minion_annotations_geneID_short.csv")
dim(annot)
names(annot)
probes = names(datExpr0)
probes2annot = match(probes, annot$transcript_id)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))

# Create the starting data frame with all the annotations, GS.weight and p.GS.weight
geneInfo0 = data.frame(transcript_id = probes,
swissprot_org = annot$swissprot_org[probes2annot],
swissprot_gene = annot$swissprot_gene[probes2annot],
swissprot_function = annot$swissprot_function[probes2annot],
interproscan = annot$interproscan[probes2annot],
TFs = annot$TFs[probes2annot],
Antismash_type = annot$Antismash_type[probes2annot],
Antismash_cluster = annot$Antismash_cluster[probes2annot],
moduleColor = modulecolours,
geneTraitSignificance,
GSPvalue)

# Order modules by their significance for weight (from smallest to largets p-val). I uses MYRO weight
modOrder = order(-abs(cor(MEs, weight, use = "p")))
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
oldNames = names(geneInfo0)
geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
MMPvalue[, modOrder[mod]]);
names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight))
geneInfo = geneInfo0[geneOrder, ]
#This data frame can be written into a text-format spreadsheet, for example by
write.csv(geneInfo, file = "geneInfo_TJ.csv")
```

# Beet_Molasses

```r
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$Beet_Molasses)
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr0, weight, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

#Summary output of network analysis results

names(datExpr0)

#names(datExpr0)[modulecolours=="green"]

annot = read.csv(file = "WT_minion_annotations_geneID_short.csv")
dim(annot)
names(annot)
probes = names(datExpr0)
probes2annot = match(probes, annot$transcript_id)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))

# Create the starting data frame with all the annotations, GS.weight and p.GS.weight
geneInfo0 = data.frame(transcript_id = probes,
swissprot_org = annot$swissprot_org[probes2annot],
swissprot_gene = annot$swissprot_gene[probes2annot],
swissprot_function = annot$swissprot_function[probes2annot],
interproscan = annot$interproscan[probes2annot],
TFs = annot$TFs[probes2annot],
Antismash_type = annot$Antismash_type[probes2annot],
Antismash_cluster = annot$Antismash_cluster[probes2annot],
moduleColor = modulecolours,
geneTraitSignificance,
GSPvalue)

# Order modules by their significance for weight (from smallest to largets p-val). I uses MYRO weight
modOrder = order(-abs(cor(MEs, weight, use = "p")))
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
oldNames = names(geneInfo0)
geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
MMPvalue[, modOrder[mod]]);
names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight))
geneInfo = geneInfo0[geneOrder, ]
#This data frame can be written into a text-format spreadsheet, for example by
write.csv(geneInfo, file = "geneInfo_MOL.csv")
```

# Maltose

```r
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$Maltose)
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr0, weight, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

#Summary output of network analysis results

names(datExpr0)

#names(datExpr0)[modulecolours=="green"]

annot = read.csv(file = "WT_minion_annotations_geneID_short.csv")
dim(annot)
names(annot)
probes = names(datExpr0)
probes2annot = match(probes, annot$transcript_id)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))

# Create the starting data frame with all the annotations, GS.weight and p.GS.weight
geneInfo0 = data.frame(transcript_id = probes,
swissprot_org = annot$swissprot_org[probes2annot],
swissprot_gene = annot$swissprot_gene[probes2annot],
swissprot_function = annot$swissprot_function[probes2annot],
interproscan = annot$interproscan[probes2annot],
TFs = annot$TFs[probes2annot],
Antismash_type = annot$Antismash_type[probes2annot],
Antismash_cluster = annot$Antismash_cluster[probes2annot],
moduleColor = modulecolours,
geneTraitSignificance,
GSPvalue)

# Order modules by their significance for weight (from smallest to largets p-val). I uses MYRO weight
modOrder = order(-abs(cor(MEs, weight, use = "p")))
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
oldNames = names(geneInfo0)
geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
MMPvalue[, modOrder[mod]]);
names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight))
geneInfo = geneInfo0[geneOrder, ]
#This data frame can be written into a text-format spreadsheet, for example by
write.csv(geneInfo, file = "geneInfo_Maltose.csv")
```


# Fructose_55percent

```r
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$Fructose_55percent)
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr0, weight, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

#Summary output of network analysis results

names(datExpr0)

#names(datExpr0)[modulecolours=="green"]

annot = read.csv(file = "WT_minion_annotations_geneID_short.csv")
dim(annot)
names(annot)
probes = names(datExpr0)
probes2annot = match(probes, annot$transcript_id)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))

# Create the starting data frame with all the annotations, GS.weight and p.GS.weight
geneInfo0 = data.frame(transcript_id = probes,
swissprot_org = annot$swissprot_org[probes2annot],
swissprot_gene = annot$swissprot_gene[probes2annot],
swissprot_function = annot$swissprot_function[probes2annot],
interproscan = annot$interproscan[probes2annot],
TFs = annot$TFs[probes2annot],
Antismash_type = annot$Antismash_type[probes2annot],
Antismash_cluster = annot$Antismash_cluster[probes2annot],
moduleColor = modulecolours,
geneTraitSignificance,
GSPvalue)

# Order modules by their significance for weight (from smallest to largets p-val). I uses MYRO weight
modOrder = order(-abs(cor(MEs, weight, use = "p")))
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
oldNames = names(geneInfo0)
geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
MMPvalue[, modOrder[mod]]);
names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight))
geneInfo = geneInfo0[geneOrder, ]
#This data frame can be written into a text-format spreadsheet, for example by
write.csv(geneInfo, file = "geneInfo_55Fructose.csv")
```

# Refined_Glucose

```r
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$Refined_Glucose)
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr0, weight, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

#Summary output of network analysis results

names(datExpr0)

#names(datExpr0)[modulecolours=="green"]

annot = read.csv(file = "WT_minion_annotations_geneID_short.csv")
dim(annot)
names(annot)
probes = names(datExpr0)
probes2annot = match(probes, annot$transcript_id)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))

# Create the starting data frame with all the annotations, GS.weight and p.GS.weight
geneInfo0 = data.frame(transcript_id = probes,
swissprot_org = annot$swissprot_org[probes2annot],
swissprot_gene = annot$swissprot_gene[probes2annot],
swissprot_function = annot$swissprot_function[probes2annot],
interproscan = annot$interproscan[probes2annot],
TFs = annot$TFs[probes2annot],
Antismash_type = annot$Antismash_type[probes2annot],
Antismash_cluster = annot$Antismash_cluster[probes2annot],
moduleColor = modulecolours,
geneTraitSignificance,
GSPvalue)

# Order modules by their significance for weight (from smallest to largets p-val). I uses MYRO weight
modOrder = order(-abs(cor(MEs, weight, use = "p")))
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
oldNames = names(geneInfo0)
geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
MMPvalue[, modOrder[mod]]);
names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight))
geneInfo = geneInfo0[geneOrder, ]
#This data frame can be written into a text-format spreadsheet, for example by
write.csv(geneInfo, file = "geneInfo_RefinedGlu.csv")
```

# Unrefined_Glucose

```r
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$Unrefined_Glucose)
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr0, weight, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

#Summary output of network analysis results

names(datExpr0)

#names(datExpr0)[modulecolours=="green"]

annot = read.csv(file = "WT_minion_annotations_geneID_short.csv")
dim(annot)
names(annot)
probes = names(datExpr0)
probes2annot = match(probes, annot$transcript_id)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))

# Create the starting data frame with all the annotations, GS.weight and p.GS.weight
geneInfo0 = data.frame(transcript_id = probes,
swissprot_org = annot$swissprot_org[probes2annot],
swissprot_gene = annot$swissprot_gene[probes2annot],
swissprot_function = annot$swissprot_function[probes2annot],
interproscan = annot$interproscan[probes2annot],
TFs = annot$TFs[probes2annot],
Antismash_type = annot$Antismash_type[probes2annot],
Antismash_cluster = annot$Antismash_cluster[probes2annot],
moduleColor = modulecolours,
geneTraitSignificance,
GSPvalue)

# Order modules by their significance for weight (from smallest to largets p-val). I uses MYRO weight
modOrder = order(-abs(cor(MEs, weight, use = "p")))
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
oldNames = names(geneInfo0)
geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
MMPvalue[, modOrder[mod]]);
names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight))
geneInfo = geneInfo0[geneOrder, ]
#This data frame can be written into a text-format spreadsheet, for example by
write.csv(geneInfo, file = "geneInfo_UnrefinedGlu.csv")
```


# Scatter plots

```r
condition="Unrefined Glucose"
module = "lightcyan"
column = match(module, modNames)
moduleGenes = modulecolours==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = paste("GS for", condition),
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "lightblue")

module = "orange"
column = match(module, modNames)
moduleGenes = modulecolours==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = paste("GS for", condition),
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "darkgrey"
column = match(module, modNames)
moduleGenes = modulecolours==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = paste("GS for", condition),
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "grey60"
column = match(module, modNames)
moduleGenes = modulecolours==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = paste("GS for", condition),
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "tan"
column = match(module, modNames)
moduleGenes = modulecolours==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = paste("GS for", condition),
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "blue"
column = match(module, modNames)
moduleGenes = modulecolours==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = paste("GS for", condition),
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "lightyellow"
column = match(module, modNames)
moduleGenes = modulecolours==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = paste("GS for", condition),
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "yellow")

module = "cyan"
column = match(module, modNames)
moduleGenes = modulecolours==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = paste("GS for", condition),
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "darkgreen"
column = match(module, modNames)
moduleGenes = modulecolours==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = paste("GS for", condition),
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "darkred"
column = match(module, modNames)
moduleGenes = modulecolours==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = paste("GS for", condition),
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
```






# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr0, power = 6)

```r

#====================================================================================
# Export data to VisANT
#====================================================================================
# Select module
module = "blue"
# Select module probes
probes = names(datExpr0)
inModule = (modulecolours==module)
modProbes = probes[inModule]
#modProbes2 = substring(modProbes,1,25)
# Select the corresponding Topological Overlap
# modTOM = tom[inModule, inModule]
# dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
nTop = 50
IMConn = softConnectivity(datExpr0[, modProbes])
top = (rank(-IMConn) <= nTop)
#vis = exportNetworkToVisANT(modTOM[top, top], file = paste("VisANTInput-", module, "-top30.txt", sep=""), weighted = TRUE, threshold = 0
 #probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol)
)
# Select top 30 interconnected genes
topGenesblue<-data.frame(IMConn,modProbes)[order(-IMConn),]
topGenesblue<-topGenesblue[1:50,]
names(topGenesblue)<-c("IMConnectivity","Genes")

# Select module
module = "tan"
# Select module probes
probes = names(datExpr0)
inModule = (modulecolours==module)
modProbes = probes[inModule]
#modProbes2 = substring(modProbes,1,25)
# Select the corresponding Topological Overlap
# modTOM = tom[inModule, inModule]
# dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
nTop = 50
IMConn = softConnectivity(datExpr0[, modProbes])
top = (rank(-IMConn) <= nTop)
#vis = exportNetworkToVisANT(modTOM[top, top], file = paste("VisANTInput-", module, "-top30.txt", sep=""), weighted = TRUE, threshold = 0
 #probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol)
# Select top 30 interconnected genes
topGenestan<-data.frame(IMConn,modProbes)[order(-IMConn),]
topGenestan<-topGenestan[1:50,]
names(topGenestan)<-c("IMConnectivity","Genes")

# Select module
module = "grey60"
# Select module probes
probes = names(datExpr0)
inModule = (modulecolours==module)
modProbes = probes[inModule]
#modProbes2 = substring(modProbes,1,25)
# Select the corresponding Topological Overlap
# modTOM = tom[inModule, inModule]
# dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
nTop = 50
IMConn = softConnectivity(datExpr0[, modProbes])
top = (rank(-IMConn) <= nTop)
#vis = exportNetworkToVisANT(modTOM[top, top], file = paste("VisANTInput-", module, "-top30.txt", sep=""), weighted = TRUE, threshold = 0
 #probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol)
# Select top 30 interconnected genes
topGenesgrey60<-data.frame(IMConn,modProbes)[order(-IMConn),]
topGenesgrey60<-topGenesgrey60[1:50,]
names(topGenesgrey60)<-c("IMConnectivity","Genes")

# Select module
module = "lightcyan"
# Select module probes
probes = names(datExpr0)
inModule = (modulecolours==module)
modProbes = probes[inModule]
#modProbes2 = substring(modProbes,1,25)
# Select the corresponding Topological Overlap
# modTOM = tom[inModule, inModule]
# dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
nTop = 50
IMConn = softConnectivity(datExpr0[, modProbes])
top = (rank(-IMConn) <= nTop)
#vis = exportNetworkToVisANT(modTOM[top, top], file = paste("VisANTInput-", module, "-top30.txt", sep=""), weighted = TRUE, threshold = 0
 #probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol)
# Select top 30 interconnected genes
topGeneslightcyan<-data.frame(IMConn,modProbes)[order(-IMConn),]
topGeneslightcyan<-topGeneslightcyan[1:50,]
names(topGeneslightcyan)<-c("IMConnectivity","Genes")

# Select module
module = "darkred"
# Select module probes
probes = names(datExpr0)
inModule = (modulecolours==module)
modProbes = probes[inModule]
#modProbes2 = substring(modProbes,1,25)
# Select the corresponding Topological Overlap
# modTOM = tom[inModule, inModule]
# dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
nTop = 50
IMConn = softConnectivity(datExpr0[, modProbes])
top = (rank(-IMConn) <= nTop)
#vis = exportNetworkToVisANT(modTOM[top, top], file = paste("VisANTInput-", module, "-top30.txt", sep=""), weighted = TRUE, threshold = 0
 #probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol)
# Select top 30 interconnected genes
topGenesdarkred<-data.frame(IMConn,modProbes)[order(-IMConn),]
topGenesdarkred<-topGenesdarkred[1:50,]
names(topGenesdarkred)<-c("IMConnectivity","Genes")

# Select module
module = "darkgreen"
# Select module probes
probes = names(datExpr0)
inModule = (modulecolours==module)
modProbes = probes[inModule]
#modProbes2 = substring(modProbes,1,25)
# Select the corresponding Topological Overlap
# modTOM = tom[inModule, inModule]
# dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
nTop = 50
IMConn = softConnectivity(datExpr0[, modProbes])
top = (rank(-IMConn) <= nTop)
#vis = exportNetworkToVisANT(modTOM[top, top], file = paste("VisANTInput-", module, "-top30.txt", sep=""), weighted = TRUE, threshold = 0
 #probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol)
# Select top 30 interconnected genes
topGenesdarkgreen<-data.frame(IMConn,modProbes)[order(-IMConn),]
topGenesdarkgreen<-topGenesdarkgreen[1:50,]
names(topGenesdarkgreen)<-c("IMConnectivity","Genes")

# Select module
module = "cyan"
# Select module probes
probes = names(datExpr0)
inModule = (modulecolours==module)
modProbes = probes[inModule]
#modProbes2 = substring(modProbes,1,25)
# Select the corresponding Topological Overlap
# modTOM = tom[inModule, inModule]
# dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
nTop = 50
IMConn = softConnectivity(datExpr0[, modProbes])
top = (rank(-IMConn) <= nTop)
#vis = exportNetworkToVisANT(modTOM[top, top], file = paste("VisANTInput-", module, "-top30.txt", sep=""), weighted = TRUE, threshold = 0
 #probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol)
# Select top 30 interconnected genes
topGenescyan<-data.frame(IMConn,modProbes)[order(-IMConn),]
topGenescyan<-topGenescyan[1:50,]
names(topGenescyan)<-c("IMConnectivity","Genes")


write.table(topGenesblue, file = "topGenesBlue.txt", sep="\t")
write.table(topGenestan, file = "topGenesTan.txt", sep="\t")
write.table(topGeneslightcyan, file = "topGenesLightCyan.txt", sep="\t")
write.table(topGenescyan, file = "topGenesCyan.txt", sep="\t")
write.table(topGenesdarkred, file = "topGenesDarkred.txt", sep="\t")
write.table(topGenesdarkgreen, file = "topGenesDarkgreen.txt", sep="\t")
write.table(topGenesgrey60, file = "topGenesGrey60.txt", sep="\t")


blue<-topGenesblue[,2]
write.table(blue,file = "blueTopID.txt", row.names = FALSE,col.names = FALSE, quote=FALSE, sep="\t")
tan<-topGenestan[,2]
write.table(tan,file = "tanTopID.txt", row.names = FALSE,col.names = FALSE, quote=FALSE, sep="\t")
lightcyan<-topGeneslightcyan[,2]
write.table(lightcyan,file = "lightcyanTopID.txt", row.names = FALSE,col.names = FALSE, quote=FALSE, sep="\t")
cyan<-topGenescyan[,2]
write.table(cyan,file = "cyanTopID.txt", row.names = FALSE,col.names = FALSE, quote=FALSE, sep="\t")
darkred<-topGenesdarkred[,2]
write.table(darkred,file = "darkredTopID.txt", row.names = FALSE,col.names = FALSE, quote=FALSE, sep="\t")
darkgreen<-topGenesdarkgreen[,2]
write.table(darkgreen,file = "darkgreenTopID.txt", row.names = FALSE,col.names = FALSE, quote=FALSE, sep="\t")
grey60<-topGenesgrey60[,2]
write.table(grey60,file = "grey60TopID.txt", row.names = FALSE,col.names = FALSE, quote=FALSE, sep="\t")


# cat interpros/all_pep.fasta | sed 's/.t1//g' > interpros/all_pep_gene.fasta

# cat interpros/all_pep_gene.fasta | grep -A 1 -wFf GeneID/Fus/KEGG/F55/darkmagenta_F55_set.txt | grep -v -- "^--$" > KEGG/green_pep.fasta

seqtk/seqtk subseq KEGG/all_pep_gene.fasta blueTopID.txt > KEGG/blueTopID.fasta
seqtk/seqtk subseq KEGG/all_pep_gene.fasta tanTopID.txt > KEGG/tanTopID.fasta
seqtk/seqtk subseq KEGG/all_pep_gene.fasta lightcyanTopID.txt > KEGG/lightcyanTopID.fasta
seqtk/seqtk subseq KEGG/all_pep_gene.fasta cyanTopID.txt > KEGG/cyanTopID.fasta
seqtk/seqtk subseq KEGG/all_pep_gene.fasta darkredTopID.txt > KEGG/darkredTopID.fasta
seqtk/seqtk subseq KEGG/all_pep_gene.fasta darkgreenTopID.txt > KEGG/darkgreenTopID.fasta
seqtk/seqtk subseq KEGG/all_pep_gene.fasta grey60TopID.txt > KEGG/grey60TopID.fasta






## Output lists of genes 

```r
library(topGO)
library(GO.db)
library(biomaRt)
library(Rgraphviz)


require(ggplot2)
library(scales)


GO_relationships <-  readMappings(file = "GO/experiment_all_gene_GO_annots_geneid.tsv")

ALL_genes=as.character(data[,1])
length(ALL_genes)
head(ALL_genes)

# All genes in modules

# module = "blue"
# probes = names(datExpr0)
# inModule = (modulecolours==module)
# modProbes = probes[inModule]

# geneList=factor(as.integer(ALL_genes %in% modProbes))
# names(geneList)= ALL_genes

seqtk/seqtk subseq KEGG/all_pep_gene.fasta GO/bluetop.txt > KEGG/bluetop.txt.fasta
seqtk/seqtk subseq KEGG/all_pep_gene.fasta GO/tantop.txt > KEGG/tantop.txt.fasta
seqtk/seqtk subseq KEGG/all_pep_gene.fasta GO/blueKOtop.txt > KEGG/blueKOtop.txt.fasta
seqtk/seqtk subseq KEGG/all_pep_gene.fasta GO/tanKOtop.txt > KEGG/tanKOtop.txt.fasta
seqtk/seqtk subseq KEGG/all_pep_gene.fasta GO/darkredtop.txt > KEGG/darkredtop.txt.fasta
seqtk/seqtk subseq KEGG/all_pep_gene.fasta GO/darkgreentop.txt > KEGG/darkgreentop.txt.fasta
seqtk/seqtk subseq KEGG/all_pep_gene.fasta GO/grey60top.txt > KEGG/grey60top.txt.fasta
seqtk/seqtk subseq KEGG/all_pep_gene.fasta GO/lightcyantop.txt > KEGG/lightcyantop.txt.fasta

# MYRO

blue4GO <- read.table("GO/bluetop.txt", header = FALSE, sep = "\t", stringsAsFactors = TRUE)
blue4GO2 <- as.character(blue4GO[,1])

geneList=factor(as.integer(ALL_genes %in% blue4GO2))
names(geneList)= ALL_genes

GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, geneSel =  blue4GO2, description = "Blue_module", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment


# allRes <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = length(allGO))
# allRes

write.table(goEnrichment, file = "GO/Hub_Blue_MF_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))


GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel =  blue4GO2, description = "Blue_module_BP", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment

write.table(goEnrichment, file = "GO/Hub_Blue_BP_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))


GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList, geneSel =  blue4GO2, description = "Blue_module", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment

write.table(goEnrichment, file = "GO/Hub_Blue_CC_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')

# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))

tan4GO <- read.table("GO/tantop.txt", header = FALSE, sep = "\t", stringsAsFactors = TRUE)
tan4GO2 <- as.character(tan4GO[,1])
geneList=factor(as.integer(ALL_genes %in% tan4GO2))
names(geneList)= ALL_genes

GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, geneSel =  tan4GO2, description = "tan_module_MF", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment


# allRes <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = length(allGO))
# allRes

write.table(goEnrichment, file = "GO/Hub_tan_MF_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))


GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel =  tan4GO2, description = "tan_module_BP", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment

write.table(goEnrichment, file = "GO/Hub_tan_BP_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))


GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList, geneSel =  tan4GO2, description = "tan_module_CC", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment

write.table(goEnrichment, file = "GO/Hub_tan_CC_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')

# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))

# MYRO_KO

blue4GO <- read.table("GO/blueKOtop.txt", header = FALSE, sep = "\t", stringsAsFactors = TRUE)
blue4GO2 <- as.character(blue4GO[,1])

geneList=factor(as.integer(ALL_genes %in% blue4GO2))
names(geneList)= ALL_genes

GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, geneSel =  blue4GO2, description = "Blue_module", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment


# allRes <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = length(allGO))
# allRes

write.table(goEnrichment, file = "GO/Hub_Blue_MF_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))


GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel =  blue4GO2, description = "Blue_module_BP", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment

write.table(goEnrichment, file = "GO/Hub_Blue_BP_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))


GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList, geneSel =  blue4GO2, description = "Blue_module", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment

write.table(goEnrichment, file = "GO/Hub_Blue_CC_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')

# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))

tan4GO <- read.table("GO/tanKOtop.txt", header = FALSE, sep = "\t", stringsAsFactors = TRUE)
tan4GO2 <- as.character(tan4GO[,1])
geneList=factor(as.integer(ALL_genes %in% tan4GO2))
names(geneList)= ALL_genes

GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, geneSel =  tan4GO2, description = "tan_module_MF", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment


# allRes <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = length(allGO))
# allRes

write.table(goEnrichment, file = "GO/Hub_tan_MF_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))


GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel =  tan4GO2, description = "tan_module_BP", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment

write.table(goEnrichment, file = "GO/Hub_tan_BP_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))


GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList, geneSel =  tan4GO2, description = "tan_module_CC", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment

write.table(goEnrichment, file = "GO/Hub_tan_CC_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')

# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))


# Maltose

darkgreen4GO <- read.table("GO/darkgreentop.txt", header = FALSE, sep = "\t", stringsAsFactors = TRUE)
darkgreen4GO2 <- as.character(darkgreen4GO[,1])

geneList=factor(as.integer(ALL_genes %in% darkgreen4GO2))
names(geneList)= ALL_genes

GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, geneSel =  darkgreen4GO2, description = "Blue_module", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment


# allRes <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = length(allGO))
# allRes

write.table(goEnrichment, file = "GO/Hub_darkgreen_MF_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))


GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel =  darkgreen4GO2, description = "Blue_module_BP", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment

write.table(goEnrichment, file = "GO/Hub_darkgreen_BP_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))


GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList, geneSel =  darkgreen4GO2, description = "Blue_module", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment

write.table(goEnrichment, file = "GO/Hub_darkgreen_CC_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')

# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))

darkred4GO <- read.table("GO/darkredtop.txt", header = FALSE, sep = "\t", stringsAsFactors = TRUE)
darkred4GO2 <- as.character(darkred4GO[,1])
geneList=factor(as.integer(ALL_genes %in% darkred4GO2))
names(geneList)= ALL_genes

GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, geneSel =  darkred4GO2, description = "MF", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment


# allRes <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = length(allGO))
# allRes

write.table(goEnrichment, file = "GO/Hub_darkred_MF_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))


GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel =  darkred4GO2, description = "BP", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment

write.table(goEnrichment, file = "GO/Hub_darkred_BP_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))


GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList, geneSel =  darkred4GO2, description = "CC", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment

write.table(goEnrichment, file = "GO/Hub_darkred_CC_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')

# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))


# Molasses

lightcyan4GO <- read.table("GO/lightcyantop.txt", header = FALSE, sep = "\t", stringsAsFactors = TRUE)
lightcyan4GO2 <- as.character(lightcyan4GO[,1])

geneList=factor(as.integer(ALL_genes %in% lightcyan4GO2))
names(geneList)= ALL_genes

GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, geneSel =  lightcyan4GO2, description = "Blue_module", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment


# allRes <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = length(allGO))
# allRes

write.table(goEnrichment, file = "GO/Hub_lightcyan_MF_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))


GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel =  lightcyan4GO2, description = "Blue_module_BP", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment

write.table(goEnrichment, file = "GO/Hub_lightcyan_BP_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))


GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList, geneSel =  lightcyan4GO2, description = "Blue_module", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment

write.table(goEnrichment, file = "GO/Hub_lightcyan_CC_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')

# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))

# F55

grey604GO <- read.table("GO/grey60top.txt", header = FALSE, sep = "\t", stringsAsFactors = TRUE)
grey604GO2 <- as.character(grey604GO[,1])
geneList=factor(as.integer(ALL_genes %in% grey604GO2))
names(geneList)= ALL_genes

GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, geneSel =  grey604GO2, description = "MF", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment


# allRes <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = length(allGO))
# allRes

write.table(goEnrichment, file = "GO/Hub_grey60_MF_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))


GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel =  grey604GO2, description = "BP", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment

write.table(goEnrichment, file = "GO/Hub_grey60_BP_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))


GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList, geneSel =  grey604GO2, description = "CC", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment

write.table(goEnrichment, file = "GO/Hub_grey60_CC_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')

# out_prefix <- paste(o, "/", "TopGO_MF", sep="")
# printGraph(GOdata, resultFisher, firstSigNodes = 10, #fn.prefix = out_prefix, 
# useInfo = "all", pdfSW = TRUE)
# length(usedGO(GOdata))







ALLMyro <- read.csv("GO/bluemyro.txt", header = TRUE, sep = "\t")
ntop <- 30
ggdata <- ALLMyro[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term))

gg2 <- ggplot(ggdata,
  aes(x = Term, y = (Expected)/(Significant), size = Significant, fill = -log10(fisher))) +

  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +

  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'GO terms in Blue Module',
    subtitle = 'Top 30 terms ordered by Fisher test p-value',
    #caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001'
    ) +
 
  #geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
    # linetype = c("dotted", "longdash", "solid"),
    # colour = c("black", "black", "black"),
    # size = c(0.5, 1.5, 3)) +

  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),

    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),

    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +

  coord_flip()

ALLMyro <- read.csv("GO/BluemyroKO.txt", header = TRUE, sep = "\t")
ntop <- 40
ggdata <- ALLMyro[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term))

gg3 <- ggplot(ggdata, aes(x = Term, y = (Expected)/(Significant), size = Significant, fill = -log10(fisher))) +

  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +

  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'GO terms in Blue Module',
    subtitle = 'Top 40 terms ordered by Fisher test p-value',
    #caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001'
    ) +
 
  #geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
    # linetype = c("dotted", "longdash", "solid"),
    # colour = c("black", "black", "black"),
    # size = c(0.5, 1.5, 3)) +

  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),

    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),

    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +

  coord_flip()




ALLMyro <- read.csv("GO/BluemyroKO.txt", header = TRUE, sep = "\t")
ntop <- 50
ggdata <- ALLMyro[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term))

gg3 <- ggplot(ggdata, aes(x = Term, y = (Expected)/(Significant), size = Significant, fill = -log10(fisher))) +

  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +

  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'GO terms in Blue Module',
    subtitle = 'Top 30 terms ordered by Fisher test p-value',
    #caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001'
    ) +
 
  #geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
    # linetype = c("dotted", "longdash", "solid"),
    # colour = c("black", "black", "black"),
    # size = c(0.5, 1.5, 3)) +

  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),

    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),

    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +

  coord_flip()

















ntop <- 600
ggdata <- goEnrichment[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term))

gg2 <- ggplot(ggdata,
  aes(x = Term, y = (Expected)/(Significant), size = Significant, fill = -log10(fisher))) +

  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +

  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'GO Molecular function',
    subtitle = 'Top 30 terms ordered by Fisher test p-value',
    #caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001'
    ) +

  #geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
    # linetype = c("dotted", "longdash", "solid"),
    # colour = c("black", "black", "black"),
    # size = c(0.5, 1.5, 3)) +

  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),

    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),

    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +

  coord_flip()

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel =  blue4GO2, description = "Blue_module", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment

GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList, geneSel =  blue4GO2, description = "Blue_module", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
resultsFisher <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
goEnrichment <- GenTable(GOdata,fisher = resultsFisher,orderBy = "fisher:ratio",topNodes = 100,numChar = 99)
goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
goEnrichment$Significant <- as.numeric(goEnrichment$Significant)
goEnrichment$Expected <- as.numeric(goEnrichment$Expected)
goEnrichment <- goEnrichment[goEnrichment$fisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","Significant","Expected","fisher")]
goEnrichment











GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel =  blue4GO2, description = "Blue_module", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
sg <- sigGenes(GOdata)
str(sg)
numSigGenes(GOdata)
resultFisher <- runTest(GOdata, algorithm="classic", statistic="fisher")
#resultFisher <- runTest(GOdata, algorithm="weight01", statistic="fisher")
resultFisher
allGO = usedGO(object = GOdata)
allRes <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = length(allGO))
allRes
write.table(allRes, file = "GO/MYRO_blue_BP_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
out_prefix <- paste(o, "/", "TopGO_BP", sep="")
printGraph(GOdata, resultFisher, firstSigNodes = 10, fn.prefix = out_prefix, useInfo = "all", pdfSW = TRUE)
length(usedGO(GOdata))


GOdata <- new("topGOdata", ontology = "CC", allGenes = gene_enrichment3, geneSel = function(p) p < 0.05, description = "Test", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
sg <- sigGenes(GOdata)
str(sg)
numSigGenes(GOdata)
#resultFisher <- runTest(GOdata, algorithm="classic", statistic="fisher")
resultFisher <- runTest(GOdata, algorithm="weight01", statistic="fisher")
resultFisher
allGO = usedGO(object = GOdata)
allRes <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = length(allGO))
allRes
write.table(allRes, file = "CC_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
out_prefix <- paste(o, "/", "TopGO_CC", sep="")
printGraph(GOdata, resultFisher, firstSigNodes = 10, fn.prefix = out_prefix, useInfo = "all", pdfSW = TRUE)
length(usedGO(GOdata))




# Extract genes from specific modules and treatmets

# Match probes in the data set to the probe IDs in the annotation file
probes = names(datExpr)

data=read.csv("iUK_data_vst_F.txt", sep="\t")
dim(data)
names(data)
# Transpose
datExpr0 = as.data.frame(t(data[, c(2:121)]))

exp_data= read.table('background_experiment1.csv', header=TRUE, sep=',')
ALL=as.character(data[,1])
ALL2=factor(as.integer(ALL))
names(geneList)= bg_genes

geneList=factor(as.integer(ALL %in% blue))
names(geneList)= ALL


ALL<-as.character(t(data[,1]))
unlist(ALL, use.names=FALSE)
probes2annot = match(probes, annot$substanceBXH)
# Get the corresponding Locuis Link IDs
allLLIDs = annot$transcript_id[probes2annot]
# $ Choose interesting modules
intModules = c("green")
for (module in intModules)
{
# Select module probes
modGenes = (modulecolours==module)
# Get their entrez ID codes
modLLIDs = allLLIDs[modGenes];
# Write them into a file
fileName = paste("MYRO", module, ".txt", sep="");
write.table(as.data.frame(modLLIDs), file = fileName,
row.names = FALSE, col.names = FALSE)
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("MYRO-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
row.names = FALSE, col.names = FALSE)
```

```bash
OutDir=./
InterProTSV=WT_minion_interproscan.tsv
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
$ProgDir/GO_table.py --interpro $InterProTSV > $OutDir/experiment_all_gene_GO_annots.tsv
cat $OutDir/experiment_all_gene_GO_annots.tsv | sed 's/.t.*//g' > $OutDir/temp1.tsv
cat $OutDir/experiment_all_gene_GO_annots.tsv | cut -f2 > $OutDir/temp2.tsv
paste $OutDir/temp1.tsv $OutDir/temp2.tsv > $OutDir/experiment_all_gene_GO_annots_geneid.tsv
cp $OutDir/experiment_all_gene_GO_annots_geneid.tsv /to/WD
rm $OutDir/temp1.tsv
rm $OutDir/temp2.tsv
```

```r

module = "cyan"
# Select module probes
probes = names(datExpr0)
inModule = (modulecolours==module)
modProbes = probes[inModule]


# Input 1
green4GO <- data.frame()
green4GO <- read.table("green4GO.txt", header = FALSE, sep = "\t", stringsAsFactors = TRUE)
gene_enrichment3 <- t(green4GO[,2])
names(gene_enrichment3) <- green4GO[,1]

green4GO <- data.frame()
green4GO <- read.table("green4GO_v2.txt", header = FALSE, sep = "\t", stringsAsFactors = TRUE)
gene_enrichment3 <- t(green4GO[,2])
names(gene_enrichment3) <- green4GO[,1]

blue4GO <- data.frame()
blue4GO <- read.table("GO/MYRO_blue.txt", header = FALSE, sep = "\t", stringsAsFactors = TRUE)
gene_enrichment <- t(blue4GO[,2])
names(gene_enrichment) <- blue4GO[,1]

# Input 2
GO_relationships <-  readMappings(file = "GO/experiment_all_gene_GO_annots_geneid.tsv")

GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, geneSel =  blue, description = "Blue_module", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata

sg <- sigGenes(GOdata)
str(sg)
numSigGenes(GOdata)

resultFisher <- runTest(GOdata, algorithm="classic", statistic="fisher")
#resultFisher <- runTest(GOdata, algorithm="weight01", statistic="fisher")
resultFisher
allGO = usedGO(object = GOdata)
allRes <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = length(allGO))
allRes
write.table(allRes, file = "GO/MYRO_blue_MF_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
out_prefix <- paste(o, "/", "TopGO_MF", sep="")
printGraph(GOdata, resultFisher, firstSigNodes = 10, fn.prefix = out_prefix, useInfo = "all", pdfSW = TRUE)
length(usedGO(GOdata))


GOdata <- new("topGOdata", ontology = "BP", allGenes = gene_enrichment, geneSel = function(p) p < 0.05, description = "Test", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
sg <- sigGenes(GOdata)
str(sg)
numSigGenes(GOdata)
resultFisher <- runTest(GOdata, algorithm="classic", statistic="fisher")
#resultFisher <- runTest(GOdata, algorithm="weight01", statistic="fisher")
resultFisher
allGO = usedGO(object = GOdata)
allRes <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = length(allGO))
allRes
write.table(allRes, file = "GO/BP_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
out_prefix <- paste(o, "/", "TopGO_BP", sep="")
printGraph(GOdata, resultFisher, firstSigNodes = 10, fn.prefix = out_prefix, useInfo = "all", pdfSW = TRUE)
length(usedGO(GOdata))


GOdata <- new("topGOdata", ontology = "CC", allGenes = gene_enrichment3, geneSel = function(p) p < 0.05, description = "Test", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
sg <- sigGenes(GOdata)
str(sg)
numSigGenes(GOdata)
#resultFisher <- runTest(GOdata, algorithm="classic", statistic="fisher")
resultFisher <- runTest(GOdata, algorithm="weight01", statistic="fisher")
resultFisher
allGO = usedGO(object = GOdata)
allRes <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = length(allGO))
allRes
write.table(allRes, file = "CC_GO_enrichment.tsv", sep = "\t")
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
out_prefix <- paste(o, "/", "TopGO_CC", sep="")
printGraph(GOdata, resultFisher, firstSigNodes = 10, fn.prefix = out_prefix, useInfo = "all", pdfSW = TRUE)
length(usedGO(GOdata))




seqtk subseq interpros/all_pep_gene.fasta MYRO/darkgreen_set.csv > MYRO/darkgreen.fasta
seqtk subseq interpros/all_pep_gene.fasta MYRO/floralwhite_set.csv > MYRO/floralwhite.fasta
seqtk subseq interpros/all_pep_gene.fasta MYRO/green_set.csv > MYRO/green.fasta

seqtk subseq interpros/all_pep_gene.fasta MYRO_KO/darkgreen_set.csv > MYRO_KO/darkgreen.fasta
seqtk subseq interpros/all_pep_gene.fasta MYRO_KO/green_set.csv > MYRO_KO/green.fasta
seqtk subseq interpros/all_pep_gene.fasta MYRO_KO/midnightblue_set.csv > MYRO_KO/midnightblue.fasta
seqtk subseq interpros/all_pep_gene.fasta MYRO_KO/yellowgreen_set.csv > MYRO_KO/yellowgreen.fasta

seqtk subseq interpros/all_pep_gene.fasta RefinedGlu/bisque_RefinedGlu_set.csv > RefinedGlu/bisque_RefinedGlu.fasta
seqtk subseq interpros/all_pep_gene.fasta RefinedGlu/brown4_RefinedGlu_set.csv > RefinedGlu/brown4_RefinedGlu.fasta
seqtk subseq interpros/all_pep_gene.fasta RefinedGlu/sienna_RefinedGlu_set.csv > RefinedGlu/sienna_RefinedGlu.fasta

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 12)
# Read in the annotation file
annot = read.csv(file = "GeneAnnotation.csv");
# Select modules
modules = c("brown", "red");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(modulecolours, module))
modProbes = probes[inModule]
modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)]
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]

dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
weighted = TRUE,
threshold = 0.02,
nodeNames = modProbes,
#altNodeNames = modGenes,
nodeAttr = modulecolours[inModule]);


transcripts <- names(datexpr0)
inmodule <- is.finite(match(modulecolours, module))
modtranscripts <- transcripts[inmodule]
modtom <- tom[inmodule, inmodule]
dimnames(modtom) <- list(modtranscripts, modtranscripts)

# Write out files for Cytoscape

edgename_start <- paste("cyt_inp_edges", module, sep = "_")
edgename <- paste(edgename_start, "txt", sep = ".")
nodename_start <- paste("cyt_inp_nodes", module, sep = "_")
nodename <- paste(nodename_start, "txt", sep = ".")

cyt <- exportNetworkToCytoscape(modtom,
edgeFile = paste(edgename, sep = "/"),
nodeFile = paste(nodename, sep = "/"),
weighted = TRUE, threshold = 0.02, nodeNames = modtranscripts,
altNodeNames = modtranscripts, nodeAttr = modulecolours[inmodule])


max <- c(4038, 1100, 312, 57, 371, 1548, 1857, 457, 49,2957,486)
> my_palette <- c("blue", "cyan", "darkgreen", "darkgrey","darkred","grey60","lightcyan","lightyellow","orange","tan","grey")
> barplot(max,col=my_palette,xlab = "Module eingene",ylab = "Number of genes")






### Maltose 

ALLMyro <- read.csv("GO/ALLmaltose.txt", header = TRUE, sep = "\t")
ntop <- 50
ggdata <- ALLMyro[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term))

gg3 <- ggplot(ggdata, aes(x = Term, y = (Expected)/(Significant), size = Significant, fill = -log10(fisher))) +

  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +

  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'GO terms in Blue Module',
    subtitle = 'Top 30 terms ordered by Fisher test p-value',
    #caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001'
    ) +
 
  #geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
    # linetype = c("dotted", "longdash", "solid"),
    # colour = c("black", "black", "black"),
    # size = c(0.5, 1.5, 3)) +

  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),

    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),

    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +

  coord_flip()