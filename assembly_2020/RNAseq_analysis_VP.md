# RNAseq data analysis

## Differential expression with DeSeq

Original data copied in GOMEZ. This was repeated to generate a vst table.

```bash
/projects/software/R-3.6.1/bin/R
```

```R
setwd("/projects/fusarium_venenatum/GOMEZ")

# Load libraries

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


# Load data from SALMON quasi mapping

# import transcript to gene mapping info
tx2gene <- read.table("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/trans2gene.txt",header=T,sep="\t")

# import quantification files

txi.reps <- tximport(paste(list.dirs("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)

# get the sample names from the folders

mysamples <- list.dirs("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2",full.names=F,recursive=F)

# summarise to gene level. This can be done in the tximport step (txOut=F), but is easier to understand in two steps.
txi.genes <- summarizeToGene(txi.reps,tx2gene)

names(txi.genes)

# set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))

# Get TPM tables
write.table(txi.genes,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/txigenes.txt",sep="\t",na="",quote=F)
write.table(txi.reps,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/txireps.txt",sep="\t",na="",quote=F)

# Read sample metadata
# Data is unordered as it is read in. This means data must be set into the same
# order as the samples were read into mysamples before integrating metadata and
# and read counts

unorderedColData <- read.table("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/FvenCarbon_RNAseq_design3.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample),])

# Group column
colData$Group <- paste0(colData$Condition,'_', colData$Timepoint)

# Group column
colData$Group2 <- paste0(colData$Condition.1,'_', colData$Timepoint)

# Define the DESeq 'GLM' model
design <- ~ Group
dds <- DESeqDataSetFromTximport(txi.genes,colData,design)

keep <- rowSums(counts(dds)) >= 50
dds <- dds[keep,]

# Library normalisation
dds <- estimateSizeFactors(dds)

# Deseq
dds <- DESeq(dds,parallel=T)

resultsNames(dds)
###
 [1] "Intercept"            "Group_C1_T1_vs_C0_T0" "Group_C1_T2_vs_C0_T0"
 [4] "Group_C1_T3_vs_C0_T0" "Group_C1_T4_vs_C0_T0" "Group_C1_T5_vs_C0_T0"
 [7] "Group_C1_T6_vs_C0_T0" "Group_C1_T7_vs_C0_T0" "Group_C2_T1_vs_C0_T0"
[10] "Group_C2_T2_vs_C0_T0" "Group_C2_T3_vs_C0_T0" "Group_C2_T4_vs_C0_T0"
[13] "Group_C2_T5_vs_C0_T0" "Group_C2_T6_vs_C0_T0" "Group_C2_T7_vs_C0_T0"
[16] "Group_C3_T1_vs_C0_T0" "Group_C3_T2_vs_C0_T0" "Group_C3_T3_vs_C0_T0"
[19] "Group_C3_T4_vs_C0_T0" "Group_C3_T5_vs_C0_T0" "Group_C3_T6_vs_C0_T0"
[22] "Group_C3_T7_vs_C0_T0" "Group_C4_T1_vs_C0_T0" "Group_C4_T2_vs_C0_T0"
[25] "Group_C4_T3_vs_C0_T0" "Group_C4_T4_vs_C0_T0" "Group_C4_T5_vs_C0_T0"
[28] "Group_C4_T6_vs_C0_T0" "Group_C4_T7_vs_C0_T0"
###


# Exploring and exporting results

res <- results(dds)
res
summary(res)

alpha <- 0.05

# res= results(dds, alpha=alpha,contrast=c("Condition","C1","C0"))
# sig.res <- subset(res,padj<=alpha)
# sig.res <- sig.res[order(sig.res$padj),]
# sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
# sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# summary(sig.res)
# # ###
# # out of 5449 with nonzero total read count
# # adjusted p-value < 0.05
# # LFC > 0 (up)       : 2551, 47%
# # LFC < 0 (down)     : 2898, 53%
# # outliers [1]       : 0, 0%
# # low counts [2]     : 0, 0%
# # (mean count < 1)
# # ###
# write.table(sig.res,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/C1_vs_C0.txt",sep="\t",na="",quote=F)
# write.table(sig.res.upregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/C1_vs_C0_up.txt",sep="\t",na="",quote=F)
# write.table(sig.res.downregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/C1_vs_C0_down.txt",sep="\t",na="",quote=F)

# res= results(dds, alpha=alpha,contrast=c("Group","C1_T2","C0_T0"))
# sig.res <- subset(res,padj<=alpha)
# sig.res <- sig.res[order(sig.res$padj),]
# sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
# sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# summary(sig.res)
# # ###
# # out of 5747 with nonzero total read count
# # adjusted p-value < 0.05
# # LFC > 0 (up)       : 2665, 46%
# # LFC < 0 (down)     : 3082, 54%
# # outliers [1]       : 0, 0%
# # low counts [2]     : 0, 0%
# # (mean count < 1)
# # ###
# write.table(sig.res,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/C2_vs_C0.txt",sep="\t",na="",quote=F)
# write.table(sig.res.upregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/C2_vs_C0_up.txt",sep="\t",na="",quote=F)
# write.table(sig.res.downregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/C2_vs_C0_down.txt",sep="\t",na="",quote=F)

# res= results(dds, alpha=alpha,contrast=c("Condition","C3","C0"))
# sig.res <- subset(res,padj<=alpha)
# sig.res <- sig.res[order(sig.res$padj),]
# sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
# sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# summary(sig.res)
# # ###
# # out of 5747 with nonzero total read count
# # adjusted p-value < 0.05
# # LFC > 0 (up)       : 2665, 46%
# # LFC < 0 (down)     : 3082, 54%
# # outliers [1]       : 0, 0%
# # low counts [2]     : 0, 0%
# # (mean count < 1)
# # ###
# write.table(sig.res,"C3_vs_C0.txt",sep="\t",na="",quote=F)
# write.table(sig.res.upregulated,"C3_vs_C0_up.txt",sep="\t",na="",quote=F)
# write.table(sig.res.downregulated,"C3_vs_C0_down.txt",sep="\t",na="",quote=F)

# # res= results(dds, alpha=alpha,contrast=c("Condition","C4","C0"))
# # sig.res <- subset(res,padj<=alpha)
# # sig.res <- sig.res[order(sig.res$padj),]
# # sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
# # sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# # summary(sig.res)
# # ###
# # out of 4192 with nonzero total read count
# # adjusted p-value < 0.05
# # LFC > 0 (up)       : 2425, 58%
# # LFC < 0 (down)     : 1767, 42%
# # outliers [1]       : 0, 0%
# # low counts [2]     : 0, 0%
# # (mean count < 0)
# # ###
# # write.table(sig.res,"C4_vs_C0.txt",sep="\t",na="",quote=F)
# # write.table(sig.res.upregulated,"C4_vs_C0_up.txt",sep="\t",na="",quote=F)
# # write.table(sig.res.downregulated,"C4_vs_C0_down.txt",sep="\t",na="",quote=F)


# Sample Distances

# These two are the same
vst1<-varianceStabilizingTransformation(dds,blind=TRUE)
write.csv(assay(vst1), file="alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/Carbon_data_vst_T.csv")
vst2<-varianceStabilizingTransformation(dds,blind=FALSE)
write.csv(assay(vst2), file="alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/Carbon_data_vst_F.csv")

pdf("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/heatmap_vst1.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst1)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst1$Condition)
colnames(sampleDistMatrix) <- paste(vst1$Timepoint)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

pdf("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/heatmap_vst2.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst2)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst2$Condition)
colnames(sampleDistMatrix) <- paste(vst2$Timepoint)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

# PCA
data <- plotPCA(vst1, intgroup=c("Timepoint","Condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Timepoint)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(vst1))) + theme_light()
coord_fixed()
ggsave("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/PCA_vst_first.pdf", pca_plot, dpi=300, height=10, width=12)

data <- plotPCA(vst2, intgroup=c("Timepoint","Condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Timepoint)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(vst2))) + theme_light()
coord_fixed()
ggsave("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/PCA_vst_second.pdf", pca_plot, dpi=300, height=10, width=12)

#Plot using rlog transformation, showing sample names:

data <- plotPCA(vst1, intgroup=c("Condition.1","Timepoint"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Timepoint, shape=Condition.1)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
geom_text_repel(aes(label=colnames(vst1))) + theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), panel.background = element_blank(),
panel.border = element_rect(colour = "black", fill = NA, size = 1),
axis.text = element_text(size = 14), axis.title = element_text(size = 18))
coord_fixed()
ggsave("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/PCA_vst_3rd.pdf", pca_plot, dpi=300, height=10, width=12)

data <- plotPCA(vst2, intgroup=c("Condition.1","Timepoint"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Timepoint, shape=Condition.1)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
geom_text_repel(aes(label=colnames(vst2))) + theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), panel.background = element_blank(),
panel.border = element_rect(colour = "black", fill = NA, size = 1),
axis.text = element_text(size = 14), axis.title = element_text(size = 18))
coord_fixed()
ggsave("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/PCA_vst_4th.pdf", pca_plot, dpi=300, height=10, width=12)
```

<!-- ```r
# Not done and not needed but I keep the commands here
# Make a table of raw counts, normalised counts and fpkm values:
raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Sample)
write.table(raw_counts,"alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Sample)
write.table(norm_counts,"alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/normalised_counts.txt",sep="\t",na="",quote=F)

# robust may be better set at false to normalise based on total counts rather than 'library normalisation factors'
fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Sample)
write.table(fpkm_counts,"alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Sample)
write.table(fpkm_counts,"alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/fpkm_counts.txt",sep="\t",na="",quote=F)


pca(experiment.table="raw_counts.txt", type="counts",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())


pca(experiment.table="normalised_counts.txt", type="counts",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())

pca(experiment.table="fpkm_counts.txt", type="FPKM",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())

pca(experiment.table="fpkm_norm_counts.txt", type="FPKM",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())


write.csv(vst, file="vst_all.csv")
write.csv(assay(vst), file="vst_all.csv") -->
£```

### Check expression of AriA and CreA

```r
# I did this for the F and T vst data.

# Rows with samples
rawdata <- read.table("Cre_AreA_T.txt",header=T,sep="\t")
# Sample, Time, OD650 columns
reshaped <- melt(rawdata, id=c("ID"), variable.name="Timepoint", value.name="vst")
write.table(reshaped, "reg_T.txt", sep="\t")

i <- read.table("reg_T.txt",header=T,sep="\t")

ggplot(i, aes(Timepoint,vst, group=ID, color=ID))+
geom_line(alpha=0.01)+
stat_summary(aes(group=ID),
fun=mean, geom="line", size=0.5)+
facet_grid(~Condition) +
xlab("Timepoints")+
ylab("vst")+
theme_bw()+
  theme(axis.text= element_text(colour="black", size=7),
        axis.title = element_text(colour = "black", size=12),
        aspect.ratio = 1, legend.title = element_blank())

i <- read.table("Area.txt",header=T,sep="\t")

ggplot(i, aes(Timepoint,vst, group=ID, color=ID))+
geom_line(alpha=0.01)+
stat_summary(aes(group=ID),
fun=mean, geom="line", size=0.5)+
facet_grid(~Condition) +
xlab("Timepoints")+
ylab("vst")+
theme_bw()+
  theme(axis.text= element_text(colour="black", size=7),
        axis.title = element_text(colour = "black", size=12),
        aspect.ratio = 1, legend.title = element_blank())

i <- read.table("Cre.txt",header=T,sep="\t")

ggplot(i, aes(Timepoint,vst, group=ID, color=ID))+
geom_line(alpha=0.01)+
stat_summary(aes(group=ID),
fun=mean, geom="line", size=0.5)+
facet_grid(~Condition) +
xlab("Timepoints")+
ylab("vst")+
theme_bw()+
  theme(axis.text= element_text(colour="black", size=7),
        axis.title = element_text(colour = "black", size=12),
        aspect.ratio = 1, legend.title = element_blank())

# Rows with samples
rawdata <- read.table("TRI_vst.txt",header=T,sep="\t")
# Sample, Time, OD650 columns
reshaped <- melt(rawdata, id=c("ID"), variable.name="Timepoint", value.name="vst")
write.table(reshaped, "TRI_res.txt", sep="\t")

i <- read.table("TRI_res.txt",header=T,sep="\t")

ggplot(i, aes(Timepoint,vst, group=ID, color=ID))+
geom_line(alpha=0.01)+
stat_summary(aes(group=ID),
fun=mean, geom="line", size=0.5)+
facet_grid(~Condition) +
xlab("Timepoints")+
ylab("vst")+
theme_bw()+
  theme(axis.text= element_text(colour="black", size=7),
        axis.title = element_text(colour = "black", size=12),
        aspect.ratio = 1, legend.title = element_blank())
```



## Interproscan

Interproscan was used to give gene models functional annotations.


```bash
# This command will split your gene fasta file and run multiple interproscan jobs.
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
for Genes in $(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta); do
echo $Genes
$ProgDir/interproscan.sh $Genes
done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following commands:

```bash
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
for Proteins in $(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta); do
Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
echo $Strain
InterProRaw=gene_pred/interproscan/F.venenatum/WT_minion/raw
$ProgDir/append_interpro.sh $Proteins $InterProRaw
done
```

## B) SwissProt

```bash
for Proteome in $(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/swissprot/$Organism/$Strain
SwissDbDir=../../dbUniprot/swissprot_2020_June
SwissDbName=uniprot_sprot
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
sbatch -p long $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
done
```

## Looking for Transcription factors

```bash
  for Interpro in $(ls gene_pred/interproscan/F.venenatum/WT_minion/*_interproscan.tsv); do
    Organism=$(echo $Interpro | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Interpro | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/transcription_factors/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
    $ProgDir/interpro2TFs.py --InterPro $Interpro > $OutDir/"$Strain"_TF_domains.tsv
    echo "total number of transcription factors"
    cat $OutDir/"$Strain"_TF_domains.tsv | cut -f1 | sort | uniq > $OutDir/"$Strain"_TF_gene_headers.txt
    cat $OutDir/"$Strain"_TF_gene_headers.txt | wc -l
    # Gene ID rather than transcript ID
    cat $OutDir/"$Strain"_TF_gene_headers.txt | sed -e "s/.t.*//g" > $OutDir/"$Strain"_TF_geneid_headers.txt
  done
```



Analysis on local machine

```R
library(GENIE3)
library(igraph)
library(RCy3)
library(Rgraphviz)

library ("reshape2")
library ("doRNG")
library ("doParallel")

source ("dynGENIE3.R")


D1<-read.table("C1_av.txt",header=T)
D2<-read.table("C2_av.txt",header=T)
D3<-read.table("C3_av.txt",header=T)
D4<-read.table("C4_av.txt",header=T)

F1 <-t(D1)
F2 <-t(D2)
F3 <-t(D3)
F4 <-t(D4)

write.table(F1, "C1_4genie.txt", sep="\t")
write.table(F2, "C2_4genie.txt", sep="\t")
write.table(F3, "C3_4genie.txt", sep="\t")
write.table(F4, "C4_4genie.txt", sep="\t")


G1<-read.table("C1_4genie.txt",header=T)
G2<-read.table("C2_4genie.txt",header=T)
G3<-read.table("C3_4genie.txt",header=T)
G4<-read.table("C4_4genie.txt",header=T)

TS1 <- read.expr.matrix("C1_av.txt",form="rows.are.genes")
TS2 <- read.expr.matrix("C2_av.txt",form="rows.are.genes")
TS3 <- read.expr.matrix("C3_av.txt",form="rows.are.genes")
TS4 <- read.expr.matrix("C4_av.txt",form="rows.are.genes")

TS1 <- read.expr.matrix("C1_4genie.txt",form="rows.are.samples")
TS2 <- read.expr.matrix("C2_4genie.txt",form="rows.are.samples")
TS3 <- read.expr.matrix("C3_4genie.txt",form="rows.are.samples")
TS4 <- read.expr.matrix("C4_4genie.txt",form="rows.are.samples")

time.points <- list(TS1[1,], TS2[1,], TS3[1,], TS4[1,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),], TS4[2:nrow(TS4),])

# Add regulators
reg<-read.table("reg.txt",header=TRUE,sep="\t")
reg2<-reg[,1]
# Regulatos final 
df <-as.character(unlist(reg2))
regulators <- df

# Run dynGENIE3
res <- dynGENIE3(TS.data,time.points, regulators=regulators)



set.seed(123)
# Use the Extra-Trees as tree-based method
tree.method <- "ET"
# Number of randomly chosen candidate regulators at each node of a tree
K <- "all"
# Number of trees per ensemble
ntrees <- 100
# Run the method with these settings
restri <- dynGENIE3(TS.data,time.points, regulators=regulators, tree.method=tree.method, K=K, ntrees=ntrees)


link.list <- get.link.list(restri$weight.matrix, threshold=0.3)
write.table(link.list, "dyGENIE3_RF_tri_01.txt", sep="\t")
# NA error on this one


# this one works. this is the previous vst data 
TS1 <- read.expr.matrix("C1_v3.txt",form="rows.are.genes")
TS2 <- read.expr.matrix("C2_v3.txt",form="rows.are.genes")
TS3 <- read.expr.matrix("C3_v3.txt",form="rows.are.genes")
TS4 <- read.expr.matrix("C4_v3.txt",form="rows.are.genes")

time.points <- list(TS1[1,], TS2[1,], TS3[1,], TS4[1,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),], TS4[2:nrow(TS4),])

# Add regulators
reg<-read.table("reg2.txt",header=TRUE,sep="\t")
reg2<-reg[,1]
# Regulatos final 
df <-as.character(unlist(reg2))
regulators <- df

# Run dynGENIE3
res <- dynGENIE3(TS.data,time.points, regulators=regulators)



set.seed(123)
# Use the Extra-Trees as tree-based method
tree.method <- "ET"
# Number of randomly chosen candidate regulators at each node of a tree
K <- "all"
# Number of trees per ensemble
ntrees <- 100
# Run the method with these settings
restri <- dynGENIE3(TS.data,time.points, regulators=regulators, tree.method=tree.method, K=K, ntrees=ntrees)

link.list <- get.link.list(restri$weight.matrix, threshold=0.5)
write.table(link.list, "dyGENIE3_RF_tri_01.txt", sep="\t")


# Add regulators
reg<-read.table("reg2.txt",header=TRUE,sep="\t")
reg2<-reg[,1]
# Regulatos final 
df <-as.character(unlist(reg2))
regulators <- df

# Run dynGENIE3
res <- dynGENIE3(TS.data,time.points, regulators=regulators)




high expressed genes 

D1<-read.table("C1_v3.txt",header=T)
D2<-read.table("C2_v3.txt",header=T)
D3<-read.table("C3_v3.txt",header=T)
D4<-read.table("C4_v3.txt",header=T)
D5<-read.table("Glucose_High_geneID.txt",header=T)
D6<-read.table("Sucrose_High_geneID.txt",header=T)
D7<-read.table("Glucose_Low_geneID.txt",header=T)
D8<-read.table("Sucrose_Low_geneID.txt",header=T)

P1<-merge(D1,D5, by.x="ID",by.y="baseMean",all.y=TRUE)
P2<-merge(D2,D6, by.x="ID",by.y="baseMean",all.y=TRUE)
P3<-merge(D3,D7, by.x="ID",by.y="baseMean",all.y=TRUE)
P4<-merge(D4,D8, by.x="ID",by.y="baseMean",all.y=TRUE)

write.table(P1, "GluHi.txt", sep="\t",quote = FALSE)
write.table(P2, "SucHi.txt", sep="\t",quote = FALSE)
write.table(P3, "GluLo.txt", sep="\t",quote = FALSE)
write.table(P4, "SucLo.txt", sep="\t",quote = FALSE)

Z1<-read.table("GluHi.txt", header=T)
Z2<-read.table("SucHi.txt", header=T)
Z3<-read.table("GluLo.txt", header=T)
Z4<-read.table("SucLo.txt", header=T)

Q1<-merge(Z1,Z2, by.x="ID",by.y="ID",all.x=TRUE, all.y=TRUE)
Q2<-merge(Z1,Z3, by.x="ID",by.y="ID",all.x=TRUE, all.y=TRUE)
Q3<-merge(Z1,Z4, by.x="ID",by.y="ID",all.x=TRUE, all.y=TRUE)

Q4<-merge(Q1,Q2, by.x="ID",by.y="ID",all.x=TRUE, all.y=TRUE)
Q5<-merge(Q3,Q4, by.x="ID",by.y="ID",all.x=TRUE, all.y=TRUE)
write.table(Q5, "GluHi_all.txt", sep="\t",quote = FALSE)

R1<-read.table("GluHi_only.txt", header=T)

P4<-merge(Z1,D8, by.x="ID",by.y="baseMean",all.y=TRUE)

# Notes. Same number of genes needed. I can do all at the same time but adding NA to the samples (not sure if this will work). otherwise, use each condition separated
# Updated. genie3 does not work with NA values. Separate conditions.

TS1 <- read.expr.matrix("GluHi_only2.txt",form="rows.are.genes")
TS2 <- read.expr.matrix("SucHi.3.txt",form="rows.are.genes")
TS3 <- read.expr.matrix("GluLo.3.txt",form="rows.are.genes")
TS4 <- read.expr.matrix("SucLo.3.txt",form="rows.are.genes")

time.points <- list(TS1[1,], TS2[1,], TS3[1,], TS4[1,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),], TS4[2:nrow(TS4),])


T12<-read.table("WT_minion_TF_gene_only_headers5.txt",header=TRUE,sep="\t")
T15<-T12[,1]
# Regulatos final 
df <-as.character(unlist(T15))
regulators <- df

set.seed(123)
# Use the Extra-Trees as tree-based method
tree.method <- "ET"
# Number of randomly chosen candidate regulators at each node of a tree
K <- "all"
# Number of trees per ensemble
ntrees <- 100
# Run the method with these settings
restri <- dynGENIE3(TS.data,time.points, regulators=regulators, tree.method=tree.method, K=K, ntrees=ntrees)

link.list <- get.link.list(restri$weight.matrix, threshold=0.5)
write.table(link.list, "dyGENIE3_RF_tri_01.txt", sep="\t")
