# RNA-Seq analysis 

## Perform qc on RNA-Seq data

```bash
# Run fastqc
    for RawData in $(ls ../../../archives/2021_camb_general/20210128_Fvenenatum_CarbonRNAseq/rawdata/*.fq.gz | grep 'FvC0'); do
        echo $RawData
        Timepoint=$(echo $RawData | rev | cut -d '/' -f1 | rev | sed -r 's/.{10}$//g' | sed -r 's/.{4}//g')
        echo $Timepoint
        Condition=$(echo $RawData | rev | cut -d '/' -f1 | rev | sed -r 's/.{12}$//g' | sed 's/Fv//g')
        echo $Condition
        OutDir=qc_rna/RNAseq/fastqc/raw/$Condition/$Timepoint
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
        sbatch -p short $ProgDir/fastqc2.sh $RawData $OutDir
    done
```

```bash
    # Run fastq-mcf
    for RNADir in $(ls -d ../../../archives/2021_camb_general/20210128_Fvenenatum_CarbonRNAseq/rawdata); do
        FileNum=$(ls $RNADir/*_1.fq.gz | grep 'FvC1T7' | wc -l)
        for num in $(seq 1 $FileNum); do
            printf "\n"
            FileF=$(ls $RNADir/*1.fq.gz | grep 'FvC1T7' | head -n $num | tail -n1)
            FileR=$(ls $RNADir/*2.fq.gz | grep 'FvC1T7' | head -n $num | tail -n1)
            echo $FileF
            echo $FileR
            #Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1.fq.gz//g')
            #echo $Sample_Name
            Timepoint=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed -r 's/.{10}$//g' | sed -r 's/.{4}//g')
            echo $Timepoint
            Condition=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed -r 's/.{12}$//g' | sed 's/Fv//g')
            echo $Condition
            OutDir=qc_rna/RNAseq/Fvenenatum_CarbonRNAseq/$Condition/$Timepoint
            echo $OutDir
            IluminaAdapters=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc/illumina_full_adapters.fa
            ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
            sbatch -p himem $ProgDir/fastq-mcf_himem.sh $FileF $FileR $IluminaAdapters RNA $OutDir
        done
    done
```

```bash
# Run fastqc
    for RawData in $(ls qc_rna/RNAseq/Fvenenatum_CarbonRNAseq/*/*/*/*.fq.gz | grep 'FvC1T7'); do
        echo $RawData
        Timepoint=$(echo $RawData | rev | cut -d '/' -f3 | rev )
        echo $Timepoint
        Condition=$(echo $RawData | rev | cut -d '/' -f4 | rev )
        echo $Condition
        OutDir=qc_rna/RNAseq/fastqc/qc_rna/$Condition/$Timepoint
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
        sbatch -p short $ProgDir/fastqc2.sh $RawData $OutDir
    done
```

## Decontamination of rRNA reads in RNAseq data

```bash
    for RNADir in $(ls -d qc_rna/RNAseq/Fvenenatum_CarbonRNAseq/C1/T7); do
        FileNum=$(ls $RNADir/F/*_1_trim.fq.gz | wc -l)
        for num in $(seq 1 $FileNum); do
            printf "\n"
            FileF=$(ls $RNADir/F/*trim.fq.gz | head -n $num | tail -n1)
            FileR=$(ls $RNADir/R/*trim.fq.gz | head -n $num | tail -n1)
            echo $FileF
            echo $FileR
            Ref=/data/scratch/gomeza/prog/bbmap/ribokmers.fa.gz
            ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
            echo $RNADir
            Strain=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
            Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
            Condition=$(echo $RNADir | rev | cut -f2 -d '/' | rev)
            echo $Condition
            Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
            echo $Sample_Name
            echo $Timepoint
            echo $Strain
            sbatch -p himem $ProgDir/bbduk.sh $Ref "$RNADir"/cleaned/$Condition/$Timepoint/$Sample_Name $FileF $FileR $ProgDir $Strain
        done
    done
```

## Salmon 


```bash
conda activate salmon
cd /projects/fusarium_venenatum

for Transcriptome in $(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.cdna.fasta); do
    Strain=$(echo $Transcriptome| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Transcriptome| rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    for RNADir in $(ls -d ../../../data/scratch/gomeza/qc_rna/RNAseq/Fvenenatum_CarbonRNAseq/C1/T7/cleaned/C*/T*/*); do
        FileNum=$(ls $RNADir/F/*_1_cleaned.fq.gz | wc -l)
        for num in $(seq 1 $FileNum); do
            printf "\n"
            FileF=$(ls $RNADir/F/*cleaned.fq.gz | head -n $num | tail -n1)
            FileR=$(ls $RNADir/R/*cleaned.fq.gz | head -n $num | tail -n1)
            echo $FileF
            echo $FileR
            Prefix=$(echo $RNADir | rev | cut -f3 -d '/' | rev)
            Timepoint=$(echo $RNADir | rev | cut -f2 -d '/' | rev)
            Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_cleaned.fq.gz//g')
            echo "$Prefix"
            echo "$Timepoint"
            echo "$Sample_Name"
            OutDir=alignment/salmon/Fvenenatum_CarbonRNAseq/$Organism/$Strain/$Prefix/$Timepoint/$Sample_Name
            ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
            sbatch -p himem $ProgDir/salmon.sh $Transcriptome $FileF $FileR $OutDir
        done
    done
done
```

Convert Salmon quasi-quanitifcations to gene counts using an awk script:

```bash
mkdir -p alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2
#This command creates a two column file with transcript_id and gene_id.
for File in $(ls alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/*/*/*/quant.sf | head -n1); do
cat $File | awk -F"\t" '{c=$1;sub(".t.*","",$1);print c,$1}' OFS="\t" > alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/trans2gene.txt
done

#This command creates a two column file with transcript_id.
#for File in $(ls alignment/salmon/*/Hg199/*/*/*/quant.sf | head -n1); do
#cat $File | awk -F"\t" '{c=$1;sub("\*","",$1);print c,$1}' OFS="\t" > alignment/salmon/N.ditissima/Hg199/DeSeq2/trans2gene3.txt
#done

# Put files in a convenient location for DeSeq.
# Analysis was not performed on Apple control samples.

for File in $(ls alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/*/*/*/quant.sf); do
Prefix=$(echo $File | cut -f8 -d '/' --output-delimiter '_')
  mkdir -p alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/$Prefix
  cp $PWD/$File alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/$Prefix/quant.sf
  # rm alignment/salmon/DeSeq2/$Prefix/quant.sf
done
```

# Differential expression with DeSeq


```bash
/projects/software/R-3.6.1/bin/R
```

```R
setwd("/projects/fusarium_venenatum")

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


# Load data from SALMON quasi mapping

# Analysis in DeSeq2 folder include all samples. Deseq_v2 does not include C0T0 samples.

# import transcript to gene mapping info
tx2gene <- read.table("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/trans2gene.txt",header=T,sep="\t")

# import quantification files
txi.reps <- tximport(paste(list.dirs("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2_v2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)

# get the sample names from the folders
mysamples <- list.dirs("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2_v2",full.names=F,recursive=F)

# summarise to gene level. This can be done in the tximport step (txOut=F), but is easier to understand in two steps.
txi.genes <- summarizeToGene(txi.reps,tx2gene)

names(txi.genes)

# set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))

# Get TPM tables
write.table(txi.genes,"alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/txigenes.txt",sep="\t",na="",quote=F)
write.table(txi.reps,"alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/txireps.txt",sep="\t",na="",quote=F)

# PCA TPMs
pca(experiment.table="alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/txigenes_TPMonly.txt", type="TPM",
legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
principal.components=c(1,2), pdf = TRUE,
output.folder=getwd())

pca(experiment.table="alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/txireps_TPMonly.txt", type="TPM",
legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
principal.components=c(1,2), pdf = TRUE,
output.folder=getwd())

# Read sample metadata
# Data is unordered as it is read in. This means data must be set into the same
# order as the samples were read into mysamples before integrating metadata and
# and read counts

unorderedColData <- read.table("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/FvenCarbon_RNAseq_design4.txt",header=T,sep="\t")
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
dds<-DESeq(dds)

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

res= results(dds, alpha=alpha,contrast=c("Condition","C1","C0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
out of 4192 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 2425, 58%
LFC < 0 (down)     : 1767, 42%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
###
write.table(sig.res,"C1_vs_C0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"C1_vs_C0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"C1_vs_C0_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Condition","C2","C0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
out of 4192 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 2425, 58%
LFC < 0 (down)     : 1767, 42%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
###
write.table(sig.res,"C2_vs_C0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"C2_vs_C0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"C2_vs_C0_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Condition","C3","C0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
out of 4192 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 2425, 58%
LFC < 0 (down)     : 1767, 42%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
###
write.table(sig.res,"C3_vs_C0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"C3_vs_C0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"C3_vs_C0_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Condition","C4","C0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
out of 4192 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 2425, 58%
LFC < 0 (down)     : 1767, 42%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
###
write.table(sig.res,"C4_vs_C0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"C4_vs_C0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"C4_vs_C0_down.txt",sep="\t",na="",quote=F)


# Sample Distances

vst1<-varianceStabilizingTransformation(dds)
vst2<-vst(dds,blind=FALSE)
vst3<-vst(dds,blind=TRUE)
pdf("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/heatmap_vst.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Condition)
colnames(sampleDistMatrix) <- paste(vst$Timepoint)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
#heatmap( sampleDistMatrix,
#  trace="none",  # turns off trace lines inside the heat map
#  col=colours, # use on color palette defined earlier
#  margins=c(12,12), # widens margins around plot
#  srtCol=45,
#  srtCol=45)
dev.off()

# Sample distances measured with rlog transformation:
 rld <- rlog(dds)
# pdf("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/heatmap_rld.pdf")
# sampleDists <- dist(t(assay(rld)))
# sampleDistMatrix <- as.matrix( sampleDists )
# rownames(sampleDistMatrix) <- paste(rld$Cultivar)
# colnames(sampleDistMatrix) <- paste(rld$Timepoint)
# colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
# dev.off()

# # MA-plot
# pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/plotMA_vst.pdf")
# plotMA(res, ylim=c(-2,2))
# dev.off()

# Plot counts
# pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/plotcounts_dds.pdf")
# plotCounts(dds, gene=which.min(res$padj), intgroup="Cultivar")
# dev.off()

# pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/plotcounts2_dds.pdf")
# plotCounts(dds, gene=which.min(res$padj), intgroup=c("Cultivar","Timepoint"))
# dev.off()

# PCA plots
pdf("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/PCA_vst_group_filtered.pdf")
plotPCA(vst,intgroup=c("Group"))
dev.off()

pdf("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2_v2/PCA_vst_FALSE.pdf")
plotPCA(vst2,intgroup=c("Group"))
dev.off()

pdf("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2_v2/PCA_vst_TRUE.pdf")
plotPCA(vst3,intgroup=c("Group"))
dev.off()

#Plot using rlog transformation:
pdf("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/PCA_rld_group.pdf")
plotPCA(rld,intgroup=c("Timepoint"))
dev.off()

#Plot using rlog transformation, showing sample names:

data <- plotPCA(vst, intgroup=c("Condition.1"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Condition.1)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(vst)))
coord_fixed()
ggsave("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/PCA_sample_names4.pdf", pca_plot, dpi=300, height=10, width=12)

data <- plotPCA(rld, intgroup=c("Condition","Timepoint"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Condition)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
coord_fixed()
ggsave("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/PCA_rld_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)


# Gene clustering plots

topVarGenes <-head(order(rowVars(assay(rld)),decreasing=TRUE),50)
pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/heatmap_new.pdf")
heatmap.2(assay(rld)[topVarGenes,],ColSideColors=c("grey","dodgerblue")[ rld$Timepoint ],scale='row',trace="none",dendrogram="column",col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
dev.off()

topVarGenes <-head(order(rowVars(assay(rld)),decreasing=TRUE),50)
pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/heatmap_new2.pdf")
colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rld$Cultivar ]
mat <- assay(rld)[ topVarGenes, ]
#mat <- mat - rowMeans(mat)
#colnames(mat) <- paste0(rld$Cultivar,"-",rld$Timepoint)
heatmap.2(mat, trace="none", col=colors, dendrogram="column",ColSideColors=sidecols,labRow=TRUE, mar=c(10,2), scale="row")
dev.off()
```

```r
# Done but not needed
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


write.csv(vst, file="vst.csv")
```

```r
# For cosistently expressed genes

topVarGenes <- head( order( rowVars( assay(vst) ), decreasing=FALSE ), 1000 ) # decreasing=TRUE for most variable genes
mat <- assay(vst)[ topVarGenes, ]
write.table(mat,"Topvar2.txt",sep="\t",na="",quote=F)

# Combine both datasets and look shared genes
T1<-read.table("Topvar1.txt",header=T,sep="\t")
T2<-read.table("Topvar2.txt",header=T,sep="\t")
T3<-merge(T1,T2, by.x="ID",by.y="ID",all.x=TRUE,all.y=TRUE) # Print all
T3<-merge(T1,T2, by.x="ID",by.y="ID",all.x=FALSE,all.y=FALSE) # Shared only
write.table(T3,"common.txt",sep="\t",na="",quote=F)
```