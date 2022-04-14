```r
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

# Read sample metadata
# Data is unordered as it is read in. This means data must be set into the same
# order as the samples were read into mysamples before integrating metadata and
# and read counts

unorderedColData <- read.table("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/FvenCarbon_RNAseq_design.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample),])

# Define the DESeq 'GLM' model
design <- ~ Condition.1
dds <- DESeqDataSetFromTximport(txi.genes,colData,design)

keep <- rowSums(counts(dds)) >= 50
dds <- dds[keep,]

# Library normalisation
dds <- estimateSizeFactors(dds)

# Deseq
dds<-DESeq(dds)

resultsNames(dds)
[1] "Intercept"                           "Condition.1_Glucose_High_vs_Control" "Condition.1_Glucose_Low_vs_Control" 
[4] "Condition.1_Sucrose_High_vs_Control" "Condition.1_Sucrose_Low_vs_Control" 

vst<-varianceStabilizingTransformation(dds,blind=FALSE)
write.csv(assay(vst), file="vst_false.csv")
vstT<-varianceStabilizingTransformation(dds,blind=TRUE)
write.csv(assay(vstT), file="vst_true.csv")


# Define the DESeq 'GLM' model
colData$Group2 <- paste0(colData$Condition.1,'_', colData$Timepoint)

design <- ~ Group2
dds <- DESeqDataSetFromTximport(txi.genes,colData,design)

keep <- rowSums(counts(dds)) >= 50
dds <- dds[keep,]

# Library normalisation
dds <- estimateSizeFactors(dds)

# Deseq
dds<-DESeq(dds)

resultsNames(dds)

 [1] "Intercept"                            "Group2_Glucose_High_T1_vs_Control_T0" "Group2_Glucose_High_T2_vs_Control_T0"
 [4] "Group2_Glucose_High_T3_vs_Control_T0" "Group2_Glucose_High_T4_vs_Control_T0" "Group2_Glucose_High_T5_vs_Control_T0"
 [7] "Group2_Glucose_High_T6_vs_Control_T0" "Group2_Glucose_High_T7_vs_Control_T0" "Group2_Glucose_Low_T1_vs_Control_T0" 
[10] "Group2_Glucose_Low_T2_vs_Control_T0"  "Group2_Glucose_Low_T3_vs_Control_T0"  "Group2_Glucose_Low_T4_vs_Control_T0" 
[13] "Group2_Glucose_Low_T5_vs_Control_T0"  "Group2_Glucose_Low_T6_vs_Control_T0"  "Group2_Glucose_Low_T7_vs_Control_T0" 
[16] "Group2_Sucrose_High_T1_vs_Control_T0" "Group2_Sucrose_High_T2_vs_Control_T0" "Group2_Sucrose_High_T3_vs_Control_T0"
[19] "Group2_Sucrose_High_T4_vs_Control_T0" "Group2_Sucrose_High_T5_vs_Control_T0" "Group2_Sucrose_High_T6_vs_Control_T0"
[22] "Group2_Sucrose_High_T7_vs_Control_T0" "Group2_Sucrose_Low_T1_vs_Control_T0"  "Group2_Sucrose_Low_T2_vs_Control_T0" 
[25] "Group2_Sucrose_Low_T3_vs_Control_T0"  "Group2_Sucrose_Low_T4_vs_Control_T0"  "Group2_Sucrose_Low_T5_vs_Control_T0" 
[28] "Group2_Sucrose_Low_T6_vs_Control_T0"  "Group2_Sucrose_Low_T7_vs_Control_T0" 

vst2F<-varianceStabilizingTransformation(dds,blind=FALSE)
write.csv(assay(vst2F), file="vst2_false.csv")
vst2T<-varianceStabilizingTransformation(dds,blind=TRUE)
write.csv(assay(vst2T), file="vst2_true.csv")





## dynGENIE3

```r
# Run on Crop diversity

# Load libraries
library(reshape2)
library(doRNG)
library(doParallel)
library(dplyr)
source ("dynGENIE3_R_C_wrapper/dynGENIE3.R")

rawF<-read.csv("vst_false.csv",header=T)

# Condition 1

# Separate replicas as independent time-series experiments
C1R1<-select(rawF,1,2,7,12,17,22,27,32,37)
C1R2<-select(rawF,1,3,8,13,18,23,28,33,38)
C1R3<-select(rawF,1,4,9,14,19,24,29,34,39)
C1R4<-select(rawF,1,5,10,15,20,25,30,35,40)
C1R5<-select(rawF,1,6,11,16,21,26,31,36,41)

#colnames(C1R1)[1]<-"gene"
colnames(C1R1)[2]<-"0"
colnames(C1R1)[3]<-"1"
colnames(C1R1)[4]<-"2"
colnames(C1R1)[5]<-"3"
colnames(C1R1)[6]<-"4"
colnames(C1R1)[7]<-"5"
colnames(C1R1)[8]<-"6"
colnames(C1R1)[9]<-"7"

#colnames(C1R2)[1]<-"gene"
colnames(C1R2)[2]<-"0"
colnames(C1R2)[3]<-"1"
colnames(C1R2)[4]<-"2"
colnames(C1R2)[5]<-"3"
colnames(C1R2)[6]<-"4"
colnames(C1R2)[7]<-"5"
colnames(C1R2)[8]<-"6"
colnames(C1R2)[9]<-"7"

colnames(C1R3)[2]<-"0"
colnames(C1R3)[3]<-"1"
colnames(C1R3)[4]<-"2"
colnames(C1R3)[5]<-"3"
colnames(C1R3)[6]<-"4"
colnames(C1R3)[7]<-"5"
colnames(C1R3)[8]<-"6"
colnames(C1R3)[9]<-"7"

colnames(C1R4)[2]<-"0"
colnames(C1R4)[3]<-"1"
colnames(C1R4)[4]<-"2"
colnames(C1R4)[5]<-"3"
colnames(C1R4)[6]<-"4"
colnames(C1R4)[7]<-"5"
colnames(C1R4)[8]<-"6"
colnames(C1R4)[9]<-"7"

colnames(C1R5)[2]<-"0"
colnames(C1R5)[3]<-"1"
colnames(C1R5)[4]<-"2"
colnames(C1R5)[5]<-"3"
colnames(C1R5)[6]<-"4"
colnames(C1R5)[7]<-"5"
colnames(C1R5)[8]<-"6"
colnames(C1R5)[9]<-"7"


# C1R1 %>% mutate_if(is.character, as.factor)
# C1R2 %>% mutate_if(is.character, as.factor)
write.table(C1R1, "C1R1.txt", sep="\t",quote = FALSE)
write.table(C1R2, "C1R2.txt", sep="\t",quote = FALSE)
write.table(C1R3, "C1R3.txt", sep="\t",quote = FALSE)
write.table(C1R4, "C1R4.txt", sep="\t",quote = FALSE)
write.table(C1R5, "C1R5.txt", sep="\t",quote = FALSE)

# Edit headers before load matrices into R with read.expr.matrix function (from dynGenie3)
TS1 <- read.expr.matrix("C1R1.txt",form="rows.are.genes")
TS2 <- read.expr.matrix("C1R2.txt",form="rows.are.genes")
TS3 <- read.expr.matrix("C1R3.txt",form="rows.are.genes")
TS4 <- read.expr.matrix("C1R4.txt",form="rows.are.genes")
TS5 <- read.expr.matrix("C1R5.txt",form="rows.are.genes")

# Time
time.points <- list(TS1[1,], TS2[1,], TS3[1,], TS4[1,], TS5[1,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),], TS4[2:nrow(TS4),], TS5[2:nrow(TS5),])

# Add regulators
TFexp<-read.table("WT_minion_TF_expression_only.txt",header=TRUE,sep="\t")
TF<-TFexp[,1]

# Default but needed
decay<-0.02
# Add to reproduce results
set.seed(123)
# Use the Extra-Trees as tree-based method
tree.method <- "ET"
# Number of randomly chosen candidate regulators at each node of a tree
K <- "all"
# Number of trees per ensemble
ntrees <- 1000
# Run the method with these settings
resall <- dynGENIE3(TS.data,time.points, regulators=TF, tree.method=tree.method, K=K, ntrees=ntrees)


# Condition 2

# Separate replicas as independent time-series experiments
C2R1<-select(rawF,1,42,47,52,57,62,67,72,77)
C2R2<-select(rawF,1,43,48,53,58,63,68,73,78)
C2R3<-select(rawF,1,44,49,54,59,64,69,74,79)
C2R4<-select(rawF,1,45,50,55,60,65,70,75,80)
C2R5<-select(rawF,1,46,51,56,61,66,71,76,81)

#colnames(C1R1)[1]<-"gene"
colnames(C2R1)[2]<-"0"
colnames(C2R1)[3]<-"1"
colnames(C2R1)[4]<-"2"
colnames(C2R1)[5]<-"3"
colnames(C2R1)[6]<-"4"
colnames(C2R1)[7]<-"5"
colnames(C2R1)[8]<-"6"
colnames(C2R1)[9]<-"7"

#colnames(C1R2)[1]<-"gene"
colnames(C2R2)[2]<-"0"
colnames(C2R2)[3]<-"1"
colnames(C2R2)[4]<-"2"
colnames(C2R2)[5]<-"3"
colnames(C2R2)[6]<-"4"
colnames(C2R2)[7]<-"5"
colnames(C2R2)[8]<-"6"
colnames(C2R2)[9]<-"7"

colnames(C2R3)[2]<-"0"
colnames(C2R3)[3]<-"1"
colnames(C2R3)[4]<-"2"
colnames(C2R3)[5]<-"3"
colnames(C2R3)[6]<-"4"
colnames(C2R3)[7]<-"5"
colnames(C2R3)[8]<-"6"
colnames(C2R3)[9]<-"7"

colnames(C2R4)[2]<-"0"
colnames(C2R4)[3]<-"1"
colnames(C2R4)[4]<-"2"
colnames(C2R4)[5]<-"3"
colnames(C2R4)[6]<-"4"
colnames(C2R4)[7]<-"5"
colnames(C2R4)[8]<-"6"
colnames(C2R4)[9]<-"7"

colnames(C2R5)[2]<-"0"
colnames(C2R5)[3]<-"1"
colnames(C2R5)[4]<-"2"
colnames(C2R5)[5]<-"3"
colnames(C2R5)[6]<-"4"
colnames(C2R5)[7]<-"5"
colnames(C2R5)[8]<-"6"
colnames(C2R5)[9]<-"7"

# C1R1 %>% mutate_if(is.character, as.factor)
# C1R2 %>% mutate_if(is.character, as.factor)
write.table(C2R1, "C2R1.txt", sep="\t",quote = FALSE)
write.table(C2R2, "C2R2.txt", sep="\t",quote = FALSE)
write.table(C2R3, "C2R3.txt", sep="\t",quote = FALSE)
write.table(C2R4, "C2R4.txt", sep="\t",quote = FALSE)
write.table(C2R5, "C2R5.txt", sep="\t",quote = FALSE)

# Edit headers before load matrices into R with read.expr.matrix function (from dynGenie3)
TS1 <- read.expr.matrix("C2R1.txt",form="rows.are.genes")
TS2 <- read.expr.matrix("C2R2.txt",form="rows.are.genes")
TS3 <- read.expr.matrix("C2R3.txt",form="rows.are.genes")
TS4 <- read.expr.matrix("C2R4.txt",form="rows.are.genes")
TS5 <- read.expr.matrix("C2R5.txt",form="rows.are.genes")

# Time
time.points <- list(TS1[1,], TS2[1,], TS3[1,], TS4[1,], TS5[1,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),], TS4[2:nrow(TS4),], TS5[2:nrow(TS5),])

# Add regulators
TFexp<-read.table("WT_minion_TF_expression_only.txt",header=TRUE,sep="\t")
TF<-TFexp[,1]

# Default but needed
decay<-0.02
# Add to reproduce results
set.seed(123)
# Use the Extra-Trees as tree-based method
tree.method <- "ET"
# Number of randomly chosen candidate regulators at each node of a tree
K <- "all"
# Number of trees per ensemble
ntrees <- 1000
# Run the method with these settings
resall <- dynGENIE3(TS.data,time.points, regulators=TF, tree.method=tree.method, K=K, ntrees=ntrees, alpha=decay)










A1<-read.table("C1_transpose.txt",header=T)
A2<-read.table("C2_transpose.txt",header=T)
A3<-read.table("C3_transpose.txt",header=T)
A4<-read.table("C4_transpose.txt",header=T)



D5<-read.table("SecmetTFs_headers.txt",header=T)
P1<-merge(D1,D5, by.x="ID",by.y="ID",all.y=TRUE)
P2<-merge(D2,D5, by.x="ID",by.y="ID",all.y=TRUE)
P3<-merge(D3,D5, by.x="ID",by.y="ID",all.y=TRUE)
P4<-merge(D4,D5, by.x="ID",by.y="ID",all.y=TRUE)
write.table(P1, "C1_v2.txt", sep="\t",quote = FALSE)
write.table(P2, "C2_v2.txt", sep="\t",quote = FALSE)
write.table(P3, "C3_v2.txt", sep="\t",quote = FALSE)
write.table(P4, "C4_v2.txt", sep="\t",quote = FALSE)
time.points <- list(TS1[1,], TS2[1,], TS3[1,], TS4[1,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),], TS4[2:nrow(TS4),])


AS1 <- read.expr.matrix("fusarium_venenatum_WD/analysis/dynGENIE3/CarbonRNAseq/v2.0/C1_transpose.txt",form="rows.are.samples")
AS2 <- read.expr.matrix("fusarium_venenatum_WD/analysis/dynGENIE3/CarbonRNAseq/v2.0/C2_transpose.txt",form="rows.are.samples")
AS3 <- read.expr.matrix("fusarium_venenatum_WD/analysis/dynGENIE3/CarbonRNAseq/v2.0/C3_transpose.txt",form="rows.are.samples")
AS4 <- read.expr.matrix("fusarium_venenatum_WD/analysis/dynGENIE3/CarbonRNAseq/v2.0/C4_transpose.txt",form="rows.are.samples")

tp <- list(AS1[1,])
AS.data <- list(AS1[2:nrow(AS1),])

time.points <- list(TS2[1,])
TS.data <- list(TS2[2:nrow(TS2),])


# Add regulators
TFexp<-read.table("WT_minion_TF_expression_only.txt",header=TRUE,sep="\t")
TF<-TFexp[,1]
# Regulatos final 
reg <-as.character(unlist(TF))
regulators <- reg

# Run dynGENIE3
res <- dynGENIE3(AS.data,tp, regulators=TF)


set.seed(123)
# Use the Extra-Trees as tree-based method
tree.method <- "ET"
# Number of randomly chosen candidate regulators at each node of a tree
K <- "all"
# Number of trees per ensemble
ntrees <- 500
# Run the method with these settings
resall <- dynGENIE3(TS.data,tp, regulators=TF, tree.method=tree.method, K=K, ntrees=ntrees)



library(reshape2)
library(doRNG)
library(doParallel)
source ("apps/prog/dynGENIE3/dynGENIE3_R_C_wrapper/dynGENIE3.R")







rsync -av *av.txt agomez@gruffalo.cropdiversity.ac.uk:/home/agomez/scratch/fusarium_venenatum_WD/analysis/ARACNe/vtest