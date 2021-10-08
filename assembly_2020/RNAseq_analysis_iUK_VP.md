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
tx2gene <- read.table("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/trans2gene.txt",header=T,sep="\t")

# import quantification files
txi.reps <- tximport(paste(list.dirs("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/temp_test", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)

# get the sample names from the folders
mysamples <- list.dirs("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/temp_test",full.names=F,recursive=F)

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
design <- ~ Media
dds <- DESeqDataSetFromTximport(txi.reps,colData,design)

# Group column
#colData$Group <- paste0(colData$Condition, colData$Sample)

# create DESeq object from count and column data
# dds <- 	DESeqDataSetFromMatrix(txi.genes,colData,~1) 
# dds <- DESeqDataSetFromTximport(txi.genes,colData,~1)	  


# add grouping factor to identify technical replicates	    
dds$groupby <- paste(dds$Media,dds$Sample,sep="_")

# sum replicates (must use same library or library size correction will go wonky)	    
dds <- collapseReplicates(dds,groupby=dds$groupby)

# normalise counts for different library size (do after collapsing replicates)
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds)) 

# define the DESeq 'GLM' model	    
# design=~Media

# add design to DESeq object	    
# design(dds) <- design # could just replace the ~1 in the first step with the design, if you really wanted to...

# Run the DESeq statistical model	    
dds <- DESeq(dds,parallel=T)


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

alpha <- 0.05

res= results(dds, alpha=alpha,contrast=c("Media","02793","02780"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/02793_vs_02780.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/02793_vs_02780_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/02793_vs_02780_down.txt",sep="\t",na="",quote=F)


res= results(dds, alpha=alpha,contrast=c("Media","10170","02780"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/10170_vs_02780.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/10170_vs_02780_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/10170_vs_02780_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Media","F55","02780"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/F55_vs_02780.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/F55_vs_02780_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/F55_vs_02780_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Media","MKO","02780"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/MKO_vs_02780.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/MKO_vs_02780_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/MKO_vs_02780_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Media","MOL","02780"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/MOL_vs_02780.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/MOL_vs_02780_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/MOL_vs_02780_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Media","MWT","02780"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/MWT_vs_02780.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/MWT_vs_02780_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/MWT_vs_02780_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Media","TJ","02780"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/TJ_vs_02780.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/TJ_vs_02780_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/TJ_vs_02780_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Media","MKO","MWT"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/MKO_vs_MWT.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/MKO_vs_MWT_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/MKO_vs_MWT_down.txt",sep="\t",na="",quote=F)

# Exploring and exporting results

# Sample Distances

vst1<-varianceStabilizingTransformation(dds,blind=TRUE)
write.csv(assay(vst1), file="alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/iUK_data_vst_T.csv")
vst2<-varianceStabilizingTransformation(dds,blind=FALSE)
write.csv(assay(vst2), file="alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/iUK_data_vst_F.csv")

pdf("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/heatmap_vst1.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst1)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst1$Media)
colnames(sampleDistMatrix) <- paste(vst1$Media)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

pdf("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/heatmap_vst2.pdf", width=12,height=12)
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
ggsave("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/PCA_vst_true.pdf", pca_plot, dpi=300, height=10, width=12)

data <- plotPCA(vst2, intgroup=c("Media"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Media)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(vst2))) + theme_light()
coord_fixed()
ggsave("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/PCA_vst_false.pdf", pca_plot, dpi=300, height=10, width=12)


rld <- rlog(dds)
#Plot using rlog transformation:
pdf("PCA_rld.pdf")
plotPCA(rld,intgroup=c("Media"))
dev.off()

# Extract genes of interest for heatmap
genes <- c("g6426.t1", "g6427.t1","g6428.t1","g6429.t1", "g6430.t1","g6431.t1", "g6432.t1","g6433.t1", "g6434.t1","g6435.t1", "g6436.t1","g6437.t1")
mat <- assay(vst2)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vst2)[c("Media")])
pheatmap(mat,color=greenred(75),annotation_col = df)

pheatmap(assay(rld)[genes,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)

heatmap.2(mat, annotation_col = anno, trace="none", key=TRUE, Colv=F, dendrogram="row", cexRow=2, cexCol=1, lwid=c(2, 6), lhei = c(2,12), density.info = "density", margins = c(1, 10))
dev.off()
```

## Analysis of DeSeq2 output

```bash
for UpFile in $(ls alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts/*_up.txt); do
DownFile=$(echo $UpFile | sed 's/_up.txt/_down.txt/g')
DegFile=$(echo $UpFile | sed 's/_up.txt/_DEGs.txt/g')
cat $UpFile $DownFile | grep -v 'baseMean' | cut -f1 | sort -u > $DegFile
echo $DegFile
cat $DegFile | wc -l
done
```

# Generating an TSV file with sequencing information

conda activate Python

```bash
for GeneGff in $(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.gff3); do
Strain=WT_minion
Organism=F.venenatum
Assembly=$(ls WT_minion_AG/repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
TFs=$(ls analysis/transcription_factors/F.venenatum/WT_minion/WT_minion_TF_domains.tsv)
InterPro=$(ls gene_pred/interproscan/F.venenatum/WT_minion/WT_minion_interproscan.tsv)
Antismash=$(ls analysis/secondary_metabolites/antismash/WT_minion_VP/WT_antismash_results_secmet_genes_corrected.tsv )
SwissProt=$(ls gene_pred/swissprot/F.venenatum/WT_minion/swissprot_vJun2020_tophit_parsed.tbl)
OutDir=analysis/annotation_tables_iUK
mkdir -p $OutDir
GeneFasta=$(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta)
Dir1=$(ls -d alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Contrasts)
DEG_Files=$(ls \
$Dir1/02793_vs_02780.txt \
$Dir1/F55_vs_02780.txt \
$Dir1/10170_vs_02780.txt \
$Dir1/MWT_vs_02780.txt \
$Dir1/MOL_vs_02780.txt \
$Dir1/MKO_vs_02780.txt \
$Dir1/TJ_vs_02780.txt \
$Dir1/MKO_vs_MWT.txt \
| sed -e "s/$/ /g" | tr -d "\n")
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Annotation_tables
python $ProgDir/build_annot_RNAseq.py  --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --DEG_files $DEG_Files --TFs $TFs --InterPro $InterPro --Antismash $Antismash --Swissprot $SwissProt > $OutDir/"$Strain"_gene_table_iUK.tsv
done
```
```bash
mv RNAseq_analysis/ RNAseq_analysis_iUK/
```

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