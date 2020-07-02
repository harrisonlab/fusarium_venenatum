

```bash
for Transcriptome in $(ls gene_pred/codingquarry_cuff_final/F.venenatum/WT_minion/final/final_genes_appended_renamed.cdna.fasta); do
Strain=$(echo $Transcriptome| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Transcriptome | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for RNADir in $(ls -d ../oldhome/groups/harrisonlab/project_files/quorn/filtered); do
FileF=$(ls $RNADir/*.1.fq | grep -e 'WTCHG_259732_224')
FileR=$(ls $RNADir/*.2.fq | grep -e 'WTCHG_259732_224')
echo $FileF
echo $FileR
Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/.1.fq//g')
echo $Sample_Name
OutDir=RNAseq_analysis/salmon_vAG/$Organism/$Strain/$Sample_Name
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
sbatch $ProgDir/salmon.sh $Transcriptome $FileF $FileR $OutDir
done
done
```
```bash


for Assembly in $(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
Query=../oldhome/groups/harrisonlab/project_files/fusarium_venenatum/analysis/genbank/g.zea_Tri5.fasta
OutDir=analysis/blast_homology/$Organism/$Strain
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
sbatch $ProgDir/blast_pipe.sh $Query protein $Assembly $OutDir
done

for BlastHits in $(ls analysis/blast_homology/*/*/*_homologs.csv); do
Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)  
Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev) 
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
HitsGff=analysis/blast_homology/$Organism/$Strain/"$Strain"_homologs.gff
Column2=BLAST_hit
NumHits=1
$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
done

SCRIPT_DIR=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
$SCRIPT_DIR/blast_self.pl ../oldhome/groups/harrisonlab/project_files/fusarium_venenatum/analysis/genbank/g.zea_Tri5.fasta blastp > g.zea_Tri5_self.csv

#-------------------------------------------------------
# 		Step 2.		simplify hits table into homolog groups
#-------------------------------------------------------

$SCRIPT_DIR/blast_parse.pl g.zea_Tri5_self.csv > g.zea_Tri5_simplified.csv

#-------------------------------------------------------
# 		Step 3.		blast queries against genome
#-------------------------------------------------------

makeblastdb -in ../oldhome/groups/harrisonlab/project_files/fusarium_venenatum/analysis/genbank/g.zea_Tri5.fasta -input_type fasta -dbtype proteins -title exons_for_blastall -parse_seqids -out exons_for_blastall

$SCRIPT_DIR/blast2csv2.pl g.zea_Tri5.fasta tblastn final_genes_appended_renamed.gene.fasta 5 > g.zea_Tri5_hits.csv

#-------------------------------------------------------
# 		Step 4.		combine the homolog group table
#					 with the blast result table
#-------------------------------------------------------

paste -d '\t' "$QUERY"_simplified.csv <(cut -f 2- "$OUTNAME"_hits.csv) > "$OUTNAME"_homologs.csv


>g6431.t1

g6427.t1;       contig_2        terpene contig_2_Cluster_10;Kind=single 
g6428.t1;       contig_2        terpene contig_2_Cluster_10;Kind=single
g6429.t1;       contig_2        terpene contig_2_Cluster_10;Kind=single
g6430.t1        contig_2        terpene contig_2_Cluster_10;Kind=single
g6431.t1;       contig_2        terpene contig_2_Cluster_10;Kind=single
g6432.t1;       contig_2        terpene contig_2_Cluster_10;Kind=single
g6433.t1        contig_2        terpene contig_2_Cluster_10;Kind=single
g6434.t1;       contig_2        terpene contig_2_Cluster_10;Kind=single
g6435.t1;       contig_2        terpene contig_2_Cluster_10;Kind=single
g6436.t1;       contig_2        terpene contig_2_Cluster_10;Kind=single



```


Convert Salmon quasi-quanitifcations to gene counts using an awk script:

```bash
mkdir -p RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2
#This command creates a two column file with transcript_id and gene_id.
for File in $(ls RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/*/quant.sf | head -n1); do
  cat $File | awk -F"\t" '{c=$1;sub(".t.*","",$1);print c,$1}' OFS="\t" > RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/trans2gene.txt
done

#This command creates a two column file with transcript_id.
#for File in $(ls alignment/salmon/*/Hg199/*/*/*/quant.sf | head -n1); do
#cat $File | awk -F"\t" '{c=$1;sub("\*","",$1);print c,$1}' OFS="\t" > alignment/salmon/N.ditissima/Hg199/DeSeq2/trans2gene3.txt
#done

# Put files in a convenient location for DeSeq.
# Analysis was not performed on Apple control samples.

for File in $(ls RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/*/quant.sf); do
  Prefix=$(echo $File | cut -f5 -d '/' --output-delimiter '_')
  mkdir -p RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/$Prefix
  cp $PWD/$File RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/$Prefix/quant.sf
  # rm alignment/salmon/DeSeq2/$Prefix/quant.sf
done
```

# Gene expression of Nd.

This analysis was done repeating the salmon alignment with the option --keepduplicates.

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

#===============================================================================
#       Load data from SALMON quasi mapping
#===============================================================================

# import transcript to gene mapping info
tx2gene <- read.table("RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/trans2gene.txt",header=T,sep="\t")

# import quantification files
txi.reps <- tximport(paste(list.dirs("RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)

# get the sample names from the folders
mysamples <- list.dirs("RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2",full.names=F,recursive=F)

# summarise to gene level. This can be done in the tximport step (txOut=F), but is easier to understand in two steps.
#txi.genes <- summarizeToGene(txi.reps,tx2gene)

# set the sample names for txi.genes
#invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))

#===============================================================================
#       Read sample metadata and annotations
#===============================================================================

# Read sample metadata
# Data is unordered as it is read in. This means data must be set into the same
# order as the samples were read into mysamples before integrating metadata and
# and read counts

unorderedColData <- read.table("RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Fven_WTminion_RNAseq_design.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample.name),])

#Add column with the media names for reference only.
colData$Media <- rep(c("02780","02793","F55","10170","MWT","MOL","MKO","TJ"),3)

# Group column
#colData$Group <- paste0(colData$Condition, colData$Sample)

# Define the DESeq 'GLM' model
design <- ~ Condition
dds <- DESeqDataSetFromTximport(txi.reps,colData,design)

# Set reference factor level, if you wish
#dds$Condition<-factor(dds$Condition, levels=c("RH1","RH2","RH3","RH4","RH5","RH6","RH7","RH8"))

# Deseq
dds<-DESeq(dds)

resultsNames(dds)
###
[1] "Intercept"            "Condition_RH2_vs_RH1" "Condition_RH3_vs_RH1"
[4] "Condition_RH4_vs_RH1" "Condition_RH5_vs_RH1" "Condition_RH6_vs_RH1"
[7] "Condition_RH7_vs_RH1" "Condition_RH8_vs_RH1"
###

#===============================================================================
#       Results
#===============================================================================

res <- results(dds)
res
summary(res)

alpha <- 0.05

res= results(dds, alpha=alpha,contrast=c("Condition","RH2","RH1"))
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
write.table(sig.res,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH2_vs_RH1.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH2_vs_RH1_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH2_vs_RH1_down.txt",sep="\t",na="",quote=F)


res= results(dds, alpha=alpha,contrast=c("Condition","RH3","RH1"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
out of 5922 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 3473, 59%
LFC < 0 (down)     : 2449, 41%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
###
write.table(sig.res,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH3_vs_RH1.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH3_vs_RH1_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH3_vs_RH1_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Condition","RH4","RH1"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
out of 7705 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 4145, 54%
LFC < 0 (down)     : 3560, 46%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
###
write.table(sig.res,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH4_vs_RH1.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH4_vs_RH1_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH4_vs_RH1_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Condition","RH5","RH1"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
out of 9062 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 4570, 50%
LFC < 0 (down)     : 4492, 50%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
###
write.table(sig.res,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH5_vs_RH1.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH5_vs_RH1_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH5_vs_RH1_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Condition","RH6","RH1"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
out of 7991 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 4029, 50%
LFC < 0 (down)     : 3962, 50%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
###
write.table(sig.res,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH6_vs_RH1.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH6_vs_RH1_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH6_vs_RH1_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Condition","RH7","RH1"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
out of 9016 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 4527, 50%
LFC < 0 (down)     : 4489, 50%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
###
write.table(sig.res,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH7_vs_RH1.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH7_vs_RH1_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH7_vs_RH1_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Condition","RH8","RH1"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
out of 7618 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 3982, 52%
LFC < 0 (down)     : 3636, 48%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
###
write.table(sig.res,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH8_vs_RH1.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH8_vs_RH1_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH8_vs_RH1_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Condition","RH8","RH1"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
out of 7618 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 3982, 52%
LFC < 0 (down)     : 3636, 48%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
###
write.table(sig.res,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH8_vs_RH1.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH8_vs_RH1_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH8_vs_RH1_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Condition","RH7","RH5"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
out of 3009 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 1456, 48%
LFC < 0 (down)     : 1553, 52%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
###
write.table(sig.res,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH7_vs_RH5.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH7_vs_RH5_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/RH7_vs_RH5_down.txt",sep="\t",na="",quote=F)

#===============================================================================
#       Exploring and exporting results
#===============================================================================

# Sample Distances

vst<-varianceStabilizingTransformation(dds)
sampleDists<-dist(t(assay(vst)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Condition)
colnames(sampleDistMatrix) <- paste(vst$Condition)
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
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Group)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

# MA-plot
pdf("MAplot.pdf")
plotMA(res, ylim=c(-2,2))
dev.off()

# Plot counts
#pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/plotcounts_dds.pdf")
plotCounts(dds, gene=which.min(res$padj), intgroup=c("Host","Sample"))
#dev.off()

# PCA plots
pdf("PCA_vst.pdf")
plotPCA(vst,intgroup=c("Condition"))
dev.off()

#Plot using rlog transformation:
pdf("PCA_rld.pdf")
plotPCA(rld,intgroup=c("Condition"))
dev.off()

#With names
pdf("PCA_rld2.pdf")
data <- plotPCA(rld, intgroup=c("Condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pca_plot<- ggplot(data, aes(PC1, PC2, color=Group, shape=Group)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
geom_text_repel(aes(label=colnames(rld))) + theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), panel.background = element_blank(),
panel.border = element_rect(colour = "black", fill = NA, size = 1),
axis.text = element_text(size = 14), axis.title = element_text(size = 18))
coord_fixed()
dev.off()

# Gene clustering plots
pdf("Top.pdf")
topVarGenes <-head(order(rowVars(assay(rld)),decreasing=TRUE),50)
heatmap.2(assay(rld)[topVarGenes,],ColSideColors=c("grey","dodgerblue")[ rld$Condition ],scale='row',trace="none",dendrogram="column",col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
dev.off()
```

## Analysis of DeSeq2 output

```bash
for UpFile in $(ls RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts/*_up.txt); do
DownFile=$(echo $UpFile | sed 's/_up.txt/_down.txt/g')
DegFile=$(echo $UpFile | sed 's/_up.txt/_DEGs.txt/g')
cat $UpFile $DownFile | grep -v 'baseMean' | cut -f1 | sort -u > $DegFile
echo $DegFile
cat $DegFile | wc -l
done
```

# Generating an TSV file with sequencing information

```bash
for GeneGff in $(ls gene_pred/codingquarry_cuff_final/F.venenatum/WT_minion/final/final_genes_appended_renamed.gff3); do
    Strain=WT_minion
    Organism=F.venenatum
    Assembly=$(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
    TFs=$(ls analysis/transcription_factors/vAG/F.venenatum/WT_minion/WT_minion_TF_domains.tsv )
    InterPro=$(ls gene_pred/interproscan/AG_medaka_47/F.venenatum/WT_minion/WT_minion_interproscan.tsv)
    #Antismash=$(ls analysis/secondary_metabolites/antismash/F.venenatum/WT/geneclusters.txt)
    Antismash=$(ls analysis/secondary_metabolites/antismash/F.venenatum/WT_minion_vAG/WT_antismash_results_secmet_genes.tsv)
    SwissProt=$(ls gene_pred/swissprot/F.venenatum/WT_minion/swissprot_vJun2020_tophit_parsed.tbl)
    OutDir=analysis/annotation_tables
    GeneFasta=$(ls gene_pred/codingquarry_cuff_final/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta)
    Dir1=$(ls -d RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/Contrasts)
    DEG_Files=$(ls \
        $Dir1/RH2_vs_RH1.txt \
        $Dir1/RH3_vs_RH1.txt \
        $Dir1/RH4_vs_RH1.txt \
        $Dir1/RH5_vs_RH1.txt \
		$Dir1/RH6_vs_RH1.txt \
		$Dir1/RH7_vs_RH1.txt \
		$Dir1/RH8_vs_RH1.txt \
        $Dir1/RH7_vs_RH5.txt \
          | sed -e "s/$/ /g" | tr -d "\n")
    ProgDir=/home/gomeza/git_repos
    $ProgDir/Nd_annotation_tables4.py  --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --DEG_files $DEG_Files --TFs $TFs --InterPro $InterPro --Antismash $Antismash --Swissprot $SwissProt > $OutDir/"$Strain"_gene_table.tsv
done


```