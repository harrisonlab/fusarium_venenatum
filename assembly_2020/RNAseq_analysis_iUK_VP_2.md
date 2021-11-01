# RNAseq iUK data analysis

## Differential expression with DeSeq

Original data copied in GOMEZ. This was repeated to generate a vst table.

```bash
/projects/software/R-3.6.1/bin/R
```

```R
setwd("/projects/fusarium_venenatum/GOMEZ")

#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library(BiocParallel)
register(MulticoreParam(12))
library(ggplot2)
library(Biostrings)
library(devtools)
library(data.table)
library(dplyr)
library(naturalsort)
library(tibble)
library(tximport)
library(rjson)
library(readr)
library(pheatmap)
library(data.table)
library(RColorBrewer)
library(gplots)
library(ggrepel)

# import transcript to gene mapping info
tx2gene <- read.table("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/data/trans2gene.txt",header=T,sep="\t")

# import quantification files
txi.reps <- tximport(paste(list.dirs("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/data", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)

# get the sample names from the folders
mysamples <- list.dirs("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/data",full.names=F,recursive=F)

# summarise to gene level. This can be done in the tximport step (txOut=F), but is easier to understand in two steps.
txi.genes <- summarizeToGene(txi.reps,tx2gene)

# set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))

# write table with TPMs
write.table(txi.genes,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/txigenes.txt",sep="\t",na="",quote=F)
write.table(txi.reps,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/txireps.txt",sep="\t",na="",quote=F)

# Read sample metadata
# Data is unordered as it is read in. This means data must be set into the same
# order as the samples were read into mysamples before integrating metadata and
# and read counts

unorderedColData <- read.table("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Fven_WTminion_RNAseq_design.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample.name),])

#Add column with the media names
colData$Media <- rep(c("02780","02793","F55","10170","MWT","MOL","MKO","TJ"),3)

# Define the DESeq 'GLM' model
design <- ~Media
dds <- DESeqDataSetFromTximport(txi.genes,colData,design)

# Group column
#colData$Group <- paste0(colData$Media,sep="_", colData$Sample)

# add grouping factor to identify technical replicates	    
dds$groupby <- paste(dds$Media,sep="_",dds$Sample)

# sum replicates (must use same library or library size correction will go wonky)	    
dds <- collapseReplicates(dds,groupby=dds$groupby)

# normalise counts for different library size (do after collapsing replicates)
#sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds)) 

# define the DESeq 'GLM' model	    
#design<-Media

# add design to DESeq object	    
#design(dds) <- design # could just replace the ~1 in the first step with the design, if you really wanted to...

# Run the DESeq statistical model	    
dds <- DESeq(dds,parallel=T)

dds2 <- collapseReplicates(dds,groupby=dds$Media)

#Pre-filtering
# This was done to extract consistently expressed genes
#keep <- rowSums(counts(dds)) >= 50
#dds <- dds[keep,]


# Set reference factor level, if you wish
#dds$Condition<-factor(dds$Condition, levels=c("RH1","RH2","RH3","RH4","RH5","RH6","RH7","RH8"))

# Deseq
# dds<-DESeq(dds)

resultsNames(dds)
###
[1] "Intercept"            "Media_02793_vs_02780" "Media_10170_vs_02780"
[4] "Media_F55_vs_02780"   "Media_MKO_vs_02780"   "Media_MOL_vs_02780"  
[7] "Media_MWT_vs_02780"   "Media_TJ_vs_02780"  
###

#===============================================================================
#       Results
#===============================================================================

res <- results(dds)
res
summary(res)

# out of 13232 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2747, 21%
# LFC < 0 (down)     : 2473, 19%

alpha <- 0.05

res= results(dds, alpha=alpha,contrast=c("Media","02793","02780"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/02793_vs_02780.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/02793_vs_02780_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/02793_vs_02780_down.txt",sep="\t",na="",quote=F)

# out of 602 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 389, 65%
# LFC < 0 (down)     : 213, 35%

res= results(dds, alpha=alpha,contrast=c("Media","10170","02780"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/10170_vs_02780.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/10170_vs_02780_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/10170_vs_02780_down.txt",sep="\t",na="",quote=F)

# out of 4374 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2453, 56%
# LFC < 0 (down)     : 1921, 44%

res= results(dds, alpha=alpha,contrast=c("Media","F55","02780"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/F55_vs_02780.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/F55_vs_02780_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/F55_vs_02780_down.txt",sep="\t",na="",quote=F)

# out of 2073 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1287, 62%
# LFC < 0 (down)     : 786, 38%

res= results(dds, alpha=alpha,contrast=c("Media","MKO","02780"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/MKO_vs_02780.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/MKO_vs_02780_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/MKO_vs_02780_down.txt",sep="\t",na="",quote=F)

# out of 7091 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3695, 52%
# LFC < 0 (down)     : 3396, 48%

res= results(dds, alpha=alpha,contrast=c("Media","MOL","02780"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/MOL_vs_02780.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/MOL_vs_02780_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/MOL_vs_02780_down.txt",sep="\t",na="",quote=F)

# out of 5145 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2721, 53%
# LFC < 0 (down)     : 2424, 47%

res= results(dds, alpha=alpha,contrast=c("Media","MWT","02780"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/MWT_vs_02780.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/MWT_vs_02780_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/MWT_vs_02780_down.txt",sep="\t",na="",quote=F)

# out of 7180 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3725, 52%
# LFC < 0 (down)     : 3455, 48%

res= results(dds, alpha=alpha,contrast=c("Media","TJ","02780"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/TJ_vs_02780.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/TJ_vs_02780_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/TJ_vs_02780_down.txt",sep="\t",na="",quote=F)

# out of 4464 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2359, 53%
# LFC < 0 (down)     : 2105, 47%

res= results(dds, alpha=alpha,contrast=c("Media","MKO","MWT"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/MKO_vs_MWT.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/MKO_vs_MWT_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/MKO_vs_MWT_down.txt",sep="\t",na="",quote=F)

# out of 140 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 52, 37%
# LFC < 0 (down)     : 88, 63%


# Exploring and exporting results

# Sample Distances

vst1<-varianceStabilizingTransformation(dds,blind=TRUE)
write.csv(assay(vst1), file="alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Tables_GeneID/iUK_data_vst_T.csv")
vst2<-varianceStabilizingTransformation(dds,blind=FALSE)
write.csv(assay(vst2), file="alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Tables_GeneID/iUK_data_vst_F.csv")

pdf("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Tables_GeneID/heatmap_vst1.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst1)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst1$Media)
colnames(sampleDistMatrix) <- paste(vst1$Media)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

pdf("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Tables_GeneID/heatmap_vst2.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst2)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst2$Media)
colnames(sampleDistMatrix) <- paste(vst2$Media)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()


# PCA
data <- plotPCA(vst1, intgroup=c("Media"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Media)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(vst1))) + theme_light()
coord_fixed()
ggsave("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Tables_GeneID/PCA_vst_true.pdf", pca_plot, dpi=300, height=10, width=12)

data <- plotPCA(vst2, intgroup=c("Media"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Media)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(vst2))) + theme_light()
coord_fixed()
ggsave("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Tables_GeneID/PCA_vst_false.pdf", pca_plot, dpi=300, height=10, width=12)


rld <- rlog(dds)
#Plot using rlog transformation:
pdf("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Tables_GeneID/PCA_rld.pdf")
plotPCA(rld,intgroup=c("Media"))
dev.off()

# Extract genes of interest for heatmap

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

vstRep2<-varianceStabilizingTransformation(dds2,blind=FALSE,fitType="local")
genes <- c("g6427","g6428","g6429", "g6430","g6431", "g6432","g6433", "g6434","g6435", "g6436")
mat <- assay(vstRep2)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vstRep2)[c("Media")])
Z<-pheatmap(mat,color=my_palette,annotation_col = anno)

vstRep1<-varianceStabilizingTransformation(dds2,blind=TRUE,fitType="local")
genes <- c("g6427","g6428","g6429", "g6430","g6431", "g6432","g6433", "g6434","g6435", "g6436")
mat <- assay(vstRep1)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vstRep1)[c("Media")])
Z<-pheatmap(mat,color=my_palette,annotation_col = anno)


vstRep2<-varianceStabilizingTransformation(dds2,blind=FALSE,fitType="parametric")
genes <- c("g6427","g6428","g6429", "g6430","g6431", "g6432","g6433", "g6434","g6435", "g6436")
mat <- assay(vstRep2)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vstRep2)[c("Media")])
pal <- wes_palette("Zissou1", 100, type = "continuous")
Z<-pheatmap(mat,color=pal,annotation_col = anno)

vstRep1<-varianceStabilizingTransformation(dds2,blind=TRUE,fitType="parametric")
genes <- c("g6427","g6428","g6429", "g6430","g6431", "g6432","g6433", "g6434","g6435", "g6436")
mat <- assay(vstRep1)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vstRep1)[c("Media")])
Z<-pheatmap(mat,color=my_palette,annotation_col = anno)

vstRep2<-varianceStabilizingTransformation(dds,blind=FALSE,fitType="local")
genes <- c("g6427","g6428","g6429", "g6430","g6431", "g6432","g6433", "g6434","g6435", "g6436")
mat <- assay(vstRep2)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vstRep2)[c("Media")])
Z<-pheatmap(mat,color=my_palette,annotation_col = anno)

vstRep1<-varianceStabilizingTransformation(dds,blind=TRUE,fitType="local")
genes <- c("g6427","g6428","g6429", "g6430","g6431", "g6432","g6433", "g6434","g6435", "g6436")
mat <- assay(vstRep1)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vstRep1)[c("Media")])
Z<-pheatmap(mat,color=my_palette,annotation_col = anno)


library(wesanderson)
gaps_row = c(5, 10, 15))
vstRep2<-varianceStabilizingTransformation(dds2,blind=FALSE,fitType="parametric")
genes <- c("g6427","g6428","g6429", "g6430","g6431", "g6432","g6433", "g6434","g6435", "g6436","g12337","g12338", "g12339","g12340", "g12341","g12342", "g12343","g12344", "g12345", "g12346","g5192","g9957","g11857"."g12859","g13028","g1506", "g1953","g2994", "g4990","g4991", "g799","g9402")
mat <- assay(vstRep2)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vstRep2)[c("Media")])
pal <- wes_palette("Zissou1", 10, type = "continuous")
Z<-pheatmap(mat,color=pal,annotation_col = anno,cutree_rows = 3, cluster_rows=FALSE, gaps_row = c(10, 20))


genes <- c("g799","g814","g1562","g1632","g1857","g1865","g1953","g2265","g2299","g2300","g2354","g2681","g2801","g3246","g3872","g4247","g4809","g5018","g5184","g5195","g5582","g6034","g6132","g6204","g6432","g6690","g6957","g7081","g7267","g7287","g8930","g9328","g9957","g11306","g11421","g11816","g11852","g11941","g11975","g12063","g12200","g13024","g13200","g13443","g13763","g13768","g13774")
mat <- assay(vstRep2)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vstRep2)[c("Media")])
pal <- wes_palette("Zissou1", 10, type = "continuous")
Z<-pheatmap(mat,color=pal,annotation_col = anno)


vstRep1<-varianceStabilizingTransformation(dds2,blind=TRUE,fitType="local")
genes <- c("g12859","g13028","g1506", "g1953","g2994", "g4990","g4991", "g799","g9402")
mat <- assay(vstRep1)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vstRep1)[c("Media")])
Z<-pheatmap(mat,color=pal,annotation_col = anno)

g12320.t1
g12337.t1
g12338.t1
g12339.t1
g12340.t1
g12341.t1
g12342.t1 - Fusarin
g12343.t1
g12344.t1
g12345.t1
g12346.t1
```

## Analysis of DeSeq2 output

```bash
for UpFile in $(ls alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/*_up.txt); do
DownFile=$(echo $UpFile | sed 's/_up.txt/_down.txt/g')
DegFile=$(echo $UpFile | sed 's/_up.txt/_DEGs.txt/g')
cat $UpFile $DownFile | grep -v 'baseMean' | cut -f1 | sort -u > $DegFile
echo $DegFile
cat $DegFile | wc -l
done
```

```
<!-- alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/02793_vs_02780_DEGs.txt
348
alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/10170_vs_02780_DEGs.txt
2517
alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/F55_vs_02780_DEGs.txt
966
alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/MKO_vs_02780_DEGs.txt
4849
alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/MKO_vs_MWT_DEGs.txt
116
alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/MOL_vs_02780_DEGs.txt
2826
alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/MWT_vs_02780_DEGs.txt
4973
alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts_GeneID/TJ_vs_02780_DEGs.txt
2303 -->
```

# Generating an TSV file with sequencing information

<!-- conda activate Python

```bash
# for GeneGff in $(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.gff3); do
# Strain=WT_minion
# Organism=F.venenatum
# Assembly=$(ls WT_minion_AG/repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
# TFs=$(ls analysis/transcription_factors/F.venenatum/WT_minion/WT_minion_TF_domains.tsv)
# InterPro=$(ls gene_pred/interproscan/F.venenatum/WT_minion/WT_minion_interproscan.tsv)
# Antismash=$(ls analysis/secondary_metabolites/antismash/WT_minion_VP/WT_antismash_results_secmet_genes_corrected.tsv )
# SwissProt=$(ls gene_pred/swissprot/F.venenatum/WT_minion/swissprot_vJun2020_tophit_parsed.tbl)
# OutDir=analysis/annotation_tables_iUK
# mkdir -p $OutDir
# GeneFasta=$(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta)
# Dir1=$(ls -d alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts)
# DEG_Files=$(ls \
# $Dir1/02793_vs_02780.txt \
# $Dir1/F55_vs_02780.txt \
# $Dir1/10170_vs_02780.txt \
# $Dir1/MWT_vs_02780.txt \
# $Dir1/MOL_vs_02780.txt \
# $Dir1/MKO_vs_02780.txt \
# $Dir1/TJ_vs_02780.txt \
# $Dir1/MKO_vs_MWT.txt \
# | sed -e "s/$/ /g" | tr -d "\n")
# ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Annotation_tables
# python $ProgDir/build_annot_RNAseq.py  --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --DEG_files $DEG_Files --TFs $TFs --InterPro $InterPro --Antismash $Antismash --Swissprot $SwissProt > $OutDir/"$Strain"_gene_table_iUK.tsv
# done
```
```bash
mv RNAseq_analysis/ RNAseq_analysis_iUK/
``` -->

# Plot genes

```bash
# Rows with samples
rawdata <- read.csv("tri5_vstF.txt",header=T,sep="\t")
# Sample, Time, OD650 columns
reshaped <- melt(rawdata, id=c("ID"), variable.name="Rep", value.name="vst")
write.table(reshaped, "TRI5res.txt", sep="\t")

Tri <- read.table("TRI5res2.txt",header=T,sep="\t")

ggplot(Tri, aes(Sample,vst, group=ID, color=ID))+
geom_line(alpha=0.01)+
stat_summary(aes(group=ID),
fun=mean, geom="line", size=0.5)+
xlab("Sample")+
ylab("vst")+
theme_bw()+
theme(axis.text= element_text(colour="black", size=7),
axis.title = element_text(colour = "black", size=12),
aspect.ratio = 1, legend.title = element_blank())

# Rows with samples
rawdata <- read.csv("tri5_vstT.txt",header=T,sep="\t")
# Sample, Time, OD650 columns
reshaped <- melt(rawdata, id=c("ID"), variable.name="Rep", value.name="vst")
write.table(reshaped, "TRI5resT.txt", sep="\t")

Tri <- read.table("TRI5resT.txt",header=T,sep="\t")

ggplot(Tri, aes(Sample,vst, group=ID, color=ID))+
geom_line(alpha=0.01)+
stat_summary(aes(group=ID),
fun=mean, geom="line", size=0.5)+
xlab("Sample")+
ylab("vst")+
theme_bw()+
theme(axis.text= element_text(colour="black", size=7),
axis.title = element_text(colour = "black", size=12),
aspect.ratio = 1, legend.title = element_blank())
```


### Coexpression analysis 


```R
# Analysis done on local machine

setwd("/projects/fusarium_venenatum/GOMEZ")

# Load libraries
library(GENIE3)
library(igraph)
library(RCy3)
library(Rgraphviz)
library(reshape2)
library(doRNG)
library(doParallel)
source ("analysis/coexpression/iUK/dynGENIE3/dynGENIE3.R")
library(CEMiTool)
library(dplyr)
```


### CEMItools


```R
# First I will run CEMiTolls with the vst and colData

# colData
# unorderedColData <- read.table("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Fven_WTminion_RNAseq_design.txt",header=T,sep="\t")
# colData <- data.frame(unorderedColData[ order(unorderedColData$Sample.name),])
# colData$Media <- rep(c("02780","02793","F55","10170","MWT","MOL","MKO","TJ"),3)

# Edit column names
colData$SampleName <- paste0(colData$Media,sep="_", colData$Sample)
colnames(colData) <- c("ID", "Condition", "Sample", "Class","Group","SampleName")

colData2<-head(colData,24)
colData3 = subset(colData2, select = -c(ID,Condition,Sample,Group))

# Datasets
t<-read.table("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Tables_GeneID/iUK_data_vst_T.txt",header=TRUE)
colnames(t) <- c("02780_A","02780_B","02780_C","02793_A","02793_B","02793_C","10170_A","10170_B","10170_C","F55_A","F55_B","F55_C","MKO_A","MKO_B","MKO_C","MOL_A","MOL_B","MOL_C","MWT_A","MWT_B","MWT_C","TJ_A","TJ_B","TJ_C")
#t2 <- t %>% mutate_all(as.numeric)
head(t[,1:4])
cem <- cemitool(t, colData3)
#"Could not specify the parameter Beta. No modules found.Unable to find parameter beta. Please check diagnostic plots with function diagnostic_report()."

# Run diagnostics to check if there are errors in the data and potential beta values. NOTE: force rewrite report
diagnostic <- diagnostic_report(cem,force=TRUE)
# Re-run with the best beta value. 
cem <- cemitool(t, colData3, set_beta=16) 

# Results
nmodules(cem) # 6 modules identified
generate_report(cem)
write_files(cem)
save_plots(cem, "all")
```

```bash
mv Plots/ analysis/coexpression/iUK/CEMiTools/vst_T_simple_gene
mv Reports/ analysis/coexpression/iUK/CEMiTools/vst_T_simple_gene
mv Tables analysis/coexpression/iUK/CEMiTools/vst_T_simple_gene
```

```r
f<-read.table("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Tables_GeneID/iUK_data_vst_F.txt",header=TRUE)
colnames(f) <- c("02780_A","02780_B","02780_C","02793_A","02793_B","02793_C","10170_A","10170_B","10170_C","F55_A","F55_B","F55_C","MKO_A","MKO_B","MKO_C","MOL_A","MOL_B","MOL_C","MWT_A","MWT_B","MWT_C","TJ_A","TJ_B","TJ_C")
#f2 <- f %>% mutate_all(as.numeric)
cem <- cemitool(f, colData3)

#"Could not specify the parameter Beta. No modules found.Unable to find parameter beta. Please check diagnostic plots with function diagnostic_report()."

# Run diagnostics to check if there are errors in the data and potential beta values. NOTE: force rewrite report
diagnostic <- diagnostic_report(cem,force=TRUE)
# Re-run with the best beta value. 
cem <- cemitool(f, colData3, set_beta=16)

# Results
nmodules(cem) # 6 modules identified
generate_report(cem)
write_files(cem)
save_plots(cem, "all")
```

```bash
mv Plots/ analysis/coexpression/iUK/CEMiTools/vst_F_simple_gene
mv Reports/ analysis/coexpression/iUK/CEMiTools/vst_F_simple_gene
mv Tables analysis/coexpression/iUK/CEMiTools/vst_F_simple_gene
```

### Genie3


```r
# Secmet antismash
# D1<-read.table("analysis/secondary_metabolites/antismash/WT_minion_VP/WT_antismash_results_secmet_genes_with_header.txt",header=T)
# less analysis/secondary_metabolites/antismash/WT_minion_VP/WT_antismash_results_secmet_genes_with_header.txt | cut -f1 | sed "s/\;//" | uniq > analysis/coexpression/iUK/Genie3/WT_minion_uniq_secmet.txt
# # All TFs no duplicate
D2<-read.table("analysis/transcription_factors/F.venenatum/WT_minion/WT_minion_TF_geneid_header_nodup.txt",header=F)
colnames(D2)[1] <- "ID"
# # Secmet and TFs
# D3<-merge(D1,D2,by.x="ID",by.y="ID",all.x=TRUE,all.y=TRUE)
# write.table(D3, "analysis/coexpression/iUK/Genie3/genes4genie3.txt", sep="\t",quote = FALSE)

# Edit: less analysis/coexpression/iUK/Genie3/genes4genie3.txt | cut -f2 | sed "s/\..*//" | uniq | sed "s/Contig/ID/" > analysis/coexpression/iUK/Genie3/genes4genie3_corrected_geneID.txt
D9<-read.table("analysis/coexpression/iUK/Genie3/genes4genie3_corrected_geneID.txt",header=TRUE)

# Expression data
write.table(t, "analysis/coexpression/iUK/Genie3/vst_T_geneid.txt", sep="\t")
write.table(f, "analysis/coexpression/iUK/Genie3/vst_F_geneid.txt", sep="\t")

# files moved to a folder called gene_id

# Add ID column manually
D4<-read.table("analysis/coexpression/iUK/Genie3/gene_id/vst_T_geneid.txt",header=T)
D5<-read.table("analysis/coexpression/iUK/Genie3/gene_id/vst_F_geneid.txt",header=T)

# Extract expression data of Secmet and TFs
D6<-merge(D4,D9, by.x="ID",by.y="ID",all.x=FALSE,all.y=FALSE)
write.table(D6, "analysis/coexpression/iUK/Genie3/gene_id/vstT4genie3.txt", sep="\t")

D7<-merge(D5,D9, by.x="ID",by.y="ID",all.x=FALSE,all.y=FALSE)
write.table(D7, "analysis/coexpression/iUK/Genie3/gene_id/vstF4genie3.txt", sep="\t")

# TFs genes only with all expression data
D10<-merge(D4,D2, by.x="ID",by.y="ID",all.x=FALSE,all.y=FALSE)
regulators<-D10[,1]

D11<-merge(D5,D2, by.x="ID",by.y="ID",all.x=FALSE,all.y=FALSE)
regulators2<-D11[,1]

# Remove ID and duplicates manually

D12<-read.table("analysis/coexpression/iUK/Genie3/gene_id/vstT4genie3.txt",header=TRUE)

D13<-read.table("analysis/coexpression/iUK/Genie3/gene_id/vstF4genie3.txt",header=TRUE)

# GENIE3

# Added for reproducibility of results
set.seed(123)

#names(D6)[1]<-""
weightMat <- GENIE3(as.matrix(D12), regulators=regulators)

#linkList <- getLinkList(weightMat)

# List for CEMiTools
linkList <- getLinkList(weightMat, threshold=0.01)
write.table(linkList, "analysis/coexpression/iUK/Genie3/gene_id/linklist_iUK_T.txt", sep="\t")

# To visalize a network
linkList <- getLinkList(weightMat, threshold=0.025)

# edge_listsi <- linkList
# # Build graph from dataframe
# Gsi <- graph.data.frame(edge_listsi,directed = T)
# # Convert graph to adjacency matrix
# Asi <- get.adjacency(Gsi,sparse = F,attr = "weight",type = "both")
# # Build adjacency graph
# g_arasi <- graph.adjacency(Asi,mode = "directed",weighted = T)
# # Create igraph
# g.cyto <- igraph.to.graphNEL(g_arasi)

# cw = createNetworkFromGraph(graph=g.cyto)

# Very different results compared to the first run

# Added for reproducibility of results
set.seed(123)

#names(D6)[1]<-""
weightMat <- GENIE3(as.matrix(D13), regulators=regulators2)

#linkList <- getLinkList(weightMat)

# List for CEMiTools
linkList <- getLinkList(weightMat, threshold=0.01)
write.table(linkList, "analysis/coexpression/iUK/Genie3/gene_id/linklist_iUK_F.txt", sep="\t")

# To visalize a network
linkList <- getLinkList(weightMat, threshold=0.025)
```


## WGCNA


```bash
OutDir=analysis/coexpression/iUK/WGCNA/
#mkdir -p $OutDir
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Co-expression_analysis
gene_table=alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Tables_GeneID/iUK_data_vst_F.txt
cat $gene_table | awk '{ print NF}'
column_start=2
column_end=25
/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/WGCNA_script.r --gene_table $gene_table --out_dir $OutDir --column_start $column_start --column_end $column_end
# Calculate softthreshold 
max_SFT=40
/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/softthreshold_power.r --out_dir $OutDir --max_SFT $max_SFT

# 4 or 9 
SFT=4
/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/create_network.r --out_dir $OutDir --sft $SFT --min_module_size 30 --merging_threshold 0.25

# Export network
/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/export_genes.r --out_dir $OutDir --unmerge Y

/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/create_network.r --out_dir analysis/coexpression/iUK/WGCNA/b9 --sft 9 --min_module_size 30 --merging_threshold 0.25

/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/create_network.r --out_dir analysis/coexpression/iUK/WGCNA/b9 --sft 16 --min_module_size 30 --merging_threshold 0.25
```

I will run WGCNA on my local machine step by step

```r
library(WGCNA)

outdir <-setwd("./")

options(stringsAsFactors = FALSE)
# Input data
data=read.csv("iUK_data_vst_F.txt", sep="\t")
dim(data)
names(data)
# Transpose
datExpr0 = as.data.frame(t(data[, c(2:25)]))
names(datExpr0) = data$ID
rownames(datExpr0) = names(data)[c(2:25)]
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

#powers <- c(c(1:40)
# Got the next message
#Warning message:
#executing %dopar% sequentially: no parallel backend registered

# Manual added seq option. No messages now
powers <- c(c(1:40),seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

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

#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html point 6 explain an alternative method to pick the softthreshold
# These plot might suggest 4 or 9 SFT values. I have 24 samples so I will use 8

adjacency <- adjacency(datExpr0, power = 8)

# Topological Overlap Matrix (TOM)

tom <- TOMsimilarity(adjacency)
disstom <- 1 - tom

# Clustering using TOM

genetree <- hclust(as.dist(disstom), method = "average")

file <- paste(outdir, "clustering_tree.pdf", sep = "/")
pdf(file, height = 9, width = 12)
sizeGrWindow(12, 9)
plot(genetree, xlab = "", sub = "", main = "Gene clustering on TOM-based
dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

# Cut clustering tree into several modules

dynamicmods <- cutreeDynamic(dendro = genetree, distM = disstom, deepSplit = 2,
pamRespectsDendro = FALSE, minClusterSize = 30)
table(dynamicmods)


# Plot modules on clustering tree, allows sanity check of min_mod_size value

file <- paste(outdir, "clustering_tree_with_modules.pdf", sep = "/")
pdf(file, height = 9, width = 12)
dynamiccolours <- labels2colors(dynamicmods)
table(dynamiccolours)
plotDendroAndColors(genetree, dynamiccolours, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colours")
dev.off()


# Merging modules with similar expression patterns

melist <- moduleEigengenes(datExpr0, colors = dynamiccolours)
mes <- melist$eigengenes
mediss <- 1 - cor(mes)
metree <- hclust(as.dist(mediss), method = "average")
file <- paste(outdir, "clustering_tree_with_merged_modules.pdf", sep = "/")
pdf(file, height = 9, width = 12)
plot(metree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h = 0.25, col = "red")
dev.off()
merge <- mergeCloseModules(datExpr0, dynamiccolours, cutHeight = 0.25,
verbose = 3)
mergedcolours <- merge$colors
mergedmes <- merge$newMEs
dev.off()

# Plot a comparison of merged and unmerged modules
file <- paste(outdir, "clustering_tree_compare_modules.pdf", sep = "/")
pdf(file, height = 9, width = 12)
plotDendroAndColors(genetree, cbind(dynamiccolours, mergedcolours),
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

```





library(qusage)
library(GSVA)
library(org.Hs.eg.db)
library(magrittr)

vst_df <- assay(vst1) %>%
  as.data.frame() %>% # Make into a data frame
  tibble::rownames_to_column("ensembl_id")