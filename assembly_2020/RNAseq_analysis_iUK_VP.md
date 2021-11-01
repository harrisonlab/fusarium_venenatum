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
txi.reps <- tximport(paste(list.dirs("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)

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
design <- ~Media
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

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

vstRep2<-varianceStabilizingTransformation(dds2,blind=F,fitType="local")
genes <- c("g6427.t1","g6428.t1","g6429.t1", "g6430.t1","g6431.t1", "g6432.t1","g6433.t1", "g6434.t1","g6435.t1", "g6436.t1")
mat <- assay(vstRep2)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vstRep2)[c("Media")])
Z<-pheatmap(mat,color=my_palette,annotation_col = anno)

vstRep1<-varianceStabilizingTransformation(dds2,blind=TRUE)
genes <- c("g6427.t1","g6428.t1","g6429.t1", "g6430.t1","g6431.t1", "g6432.t1","g6433.t1", "g6434.t1","g6435.t1", "g6436.t1")
mat <- assay(vstRep1)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vstRep1)[c("Media")])
Z<-pheatmap(mat,color=my_palette,annotation_col = anno)


vstRep2<-varianceStabilizingTransformation(dds2,blind=FALSE,fitType="parametric")
genes <- c("g6427.t1","g6428.t1","g6429.t1", "g6430.t1","g6431.t1", "g6432.t1","g6433.t1", "g6434.t1","g6435.t1", "g6436.t1")
mat <- assay(vstRep2)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vstRep2)[c("Media")])
Z<-pheatmap(mat,color=my_palette,annotation_col = anno)

vstRep1<-varianceStabilizingTransformation(dds2,blind=TRUE,fitType="parametric")
genes <- c("g6427.t1","g6428.t1","g6429.t1", "g6430.t1","g6431.t1", "g6432.t1","g6433.t1", "g6434.t1","g6435.t1", "g6436.t1")
mat <- assay(vstRep1)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vstRep1)[c("Media")])
Z<-pheatmap(mat,color=my_palette,annotation_col = anno)
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
unorderedColData <- read.table("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Fven_WTminion_RNAseq_design.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample.name),])
colData$Media <- rep(c("02780","02793","F55","10170","MWT","MOL","MKO","TJ"),3)

# Edit column names
colData$SampleName <- paste0(colData$Media,sep="_", colData$Sample)
colnames(colData) <- c("ID", "Condition", "Sample", "Class","SampleName")

colData2<-head(colData,24)
colData3 = subset(colData2, select = -c(ID,Condition,Sample))

# Datasets
t<-read.table("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/iUK_data_vst_T.txt",header=TRUE)
colnames(t) <- c("02780_A","02780_B","02780_C","02793_A","02793_B","02793_C","F55_A","F55_B","F55_C","10170_A","10170_B","10170_C","MWT_A","MWT_B","MWT_C","MOL_A","MOL_B","MOL_C","MKO_A","MKO_B","MKO_C","TJ_A","TJ_B","TJ_C")
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
mv Plots/ analysis/coexpression/iUK/CEMiTools/vst_T_simple
mv Reports/ analysis/coexpression/iUK/CEMiTools/vst_T_simple
mv Tables analysis/coexpression/iUK/CEMiTools/vst_T_simple
```

```r
f<-read.table("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/iUK_data_vst_F.txt",header=TRUE)
colnames(f) <- c("02780_A","02780_B","02780_C","02793_A","02793_B","02793_C","F55_A","F55_B","F55_C","10170_A","10170_B","10170_C","MWT_A","MWT_B","MWT_C","MOL_A","MOL_B","MOL_C","MKO_A","MKO_B","MKO_C","TJ_A","TJ_B","TJ_C")
#f2 <- f %>% mutate_all(as.numeric)
cem <- cemitool(f, colData3)

#"Could not specify the parameter Beta. No modules found.Unable to find parameter beta. Please check diagnostic plots with function diagnostic_report()."

# Run diagnostics to check if there are errors in the data and potential beta values. NOTE: force rewrite report
diagnostic <- diagnostic_report(cem,force=TRUE)
# Re-run with the best beta value. 
cem <- cemitool(f, colData3, set_beta=16)

# Results
nmodules(cem) # 5 modules identified
generate_report(cem)
write_files(cem)
save_plots(cem, "all")
```

```bash
mv Plots/ analysis/coexpression/iUK/CEMiTools/vst_F_simple/
mv Reports/ analysis/coexpression/iUK/CEMiTools/vst_F_simple/
mv Tables analysis/coexpression/iUK/CEMiTools/vst_F_simple/
```

### Genie3


```r
# Secmet antismash
D1<-read.table("analysis/secondary_metabolites/antismash/WT_minion_VP/WT_antismash_results_secmet_genes_with_header.txt",header=T)
less analysis/secondary_metabolites/antismash/WT_minion_VP/WT_antismash_results_secmet_genes_with_header.txt | cut -f1 | sed "s/\;//" | uniq > analysis/coexpression/iUK/Genie3/WT_minion_uniq_secmet.txt
# All TFs no duplicate
D2<-read.table("analysis/transcription_factors/F.venenatum/WT_minion/WT_minion_TF_gene_headers.txt",header=F)
colnames(D2)[1] <- "ID"
# Secmet and TFs
D3<-merge(D1,D2,by.x="ID",by.y="ID",all.x=TRUE,all.y=TRUE)
write.table(D3, "analysis/coexpression/iUK/Genie3/genes4genie3.txt", sep="\t",quote = FALSE)

# Edit: less analysis/coexpression/iUK/Genie3/genes4genie3.txt | cut -f2 | sed "s/\;//" | uniq | sed "s/Contig/ID/" > analysis/coexpression/iUK/Genie3/genes4genie3_corrected.txt
D9<-read.table("analysis/coexpression/iUK/Genie3/genes4genie3_corrected.txt",header=TRUE)

# Expression data
write.table(t, "analysis/coexpression/iUK/Genie3/vst_T.txt", sep="\t")
write.table(f, "analysis/coexpression/iUK/Genie3/vst_F.txt", sep="\t")
# Add ID column manually
D4<-read.table("analysis/coexpression/iUK/Genie3/vst_T.txt",header=T)
D5<-read.table("analysis/coexpression/iUK/Genie3/vst_F.txt",header=T)

# Extract expression data of Secmet and TFs
D6<-merge(D4,D9, by.x="ID",by.y="ID",all.x=FALSE,all.y=FALSE)
write.table(D6, "analysis/coexpression/iUK/Genie3/vstT4genie3.txt", sep="\t")

D6<-merge(D5,D9, by.x="ID",by.y="ID",all.x=FALSE,all.y=FALSE)
write.table(D6, "analysis/coexpression/iUK/Genie3/vstF4genie3.txt", sep="\t")

# TFs genes only with all expression data
D7<-merge(D4,D2, by.x="ID",by.y="ID",all.x=FALSE,all.y=FALSE)
regulators<-D7[,1]

D7<-merge(D5,D2, by.x="ID",by.y="ID",all.x=FALSE,all.y=FALSE)
regulators<-D7[,1]

# Remove ID and duplicates manually
D8<-read.table("analysis/coexpression/iUK/Genie3/vstT4genie3.txt",header=TRUE)

D8<-read.table("analysis/coexpression/iUK/Genie3/vstF4genie3.txt",header=TRUE)

# GENIE3

# Added for reproducibility of results
set.seed(123)

#names(D6)[1]<-""
weightMat <- GENIE3(as.matrix(D8), regulators=regulators)

#linkList <- getLinkList(weightMat)

# List for CEMiTools
linkList <- getLinkList(weightMat, threshold=0.01)
write.table(linkList, "analysis/coexpression/iUK/Genie3/linklist_iUK_F.txt", sep="\t")

# To visalize a network
linkList <- getLinkList(weightMat, threshold=0.025)

edge_listsi <- linkList
# Build graph from dataframe
Gsi <- graph.data.frame(edge_listsi,directed = T)
# Convert graph to adjacency matrix
Asi <- get.adjacency(Gsi,sparse = F,attr = "weight",type = "both")
# Build adjacency graph
g_arasi <- graph.adjacency(Asi,mode = "directed",weighted = T)
# Create igraph
g.cyto <- igraph.to.graphNEL(g_arasi)

cw = createNetworkFromGraph(graph=g.cyto)

# Very different results compared to the first run
```











TS1 <- as.data.frame.matrix("vst1_corrected.txt")
reg<-read.table("Carbon_data_vst_F.txt",header=TRUE)
reg<-read.table("vst1_corrected.txt",header=TRUE)
View(reg)
cem <- cemitool(reg)


vstT<-read.table("Carbon_data_vst_T.csv",header=TRUE,sep="\t")
cem <- cemitool(vstT, colData)

nmodules(cem)
generate_report(cem)
write_files(cem)
save_plots(cem, "all")

raw$Group <- paste0(raw$Class,'_', raw$Timepoint)


colData<-read.table("colData.txt",header=T,sep="\t")


vst_F does not work!!!

vst_T worked

I tried group and condition. both gave a thershold in 12 and not unique module in sucrose high



cem <- cemitool(vstT, colData,filter_pval=0.2) # it takes more time. One module more. No tri5


este no funciono en iUK. debo repetirlo cuando pruebe el que si funciono


reg<-read.table("vst_true4cemi.txt",header=TRUE)
reg2 <- reg %>% mutate_all(as.numeric)

reg3<-read.table("vst_false4cemi.txt",header=TRUE)
reg4 <- reg3 %>% mutate_all(as.numeric)

coldata<-read.table("colData.txt",header=T,sep="\t")

cem <- cemitool(reg, coldata)