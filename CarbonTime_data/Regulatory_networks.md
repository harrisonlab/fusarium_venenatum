# Dynamic regulatory networks

### Input files needed 

```bash
# Rename TF files. FYI, original files created in Genie3 script
cd analysis/transcription_factors
# TFs transcript id
mv WT_minion_TF_gene_headers.txt WT_minion_TF_transcript_headers.txt
# TFs gene id
mv WT_minion_TF_gene_only_headers.txt WT_minion_TF_gene_headers_ID.txt
# All TFs no duplicate no ID
mv WT_minion_TF_gene_only_headers3.txt WT_minion_TF_gene_headers_nodup_noID.txt
# All TFs no duplicate
mv WT_minion_TF_gene_only_headers4.txt WT_minion_TF_gene_headers_nodup.txt
# TFs with expression data only
mv WT_minion_TF_gene_only_headers5.txt WT_minion_TF_expression_only.txt 

conda activate perly_env
```

## dynGENIE3

```
Notes for dynGENIE3
dynGENIE3 input list of matrices with one matrix per experiment (pseudoreps). 
Create one network per condition to compare with other methods.
I recommend vst transformed data using blind=FALSE to account for differences in expression by terms in the design (e.g. media, cultivar....)
I am not using the control samples for this analysis
```

```bash
srun --partition long --time 0-02:00:00 --mem-per-cpu 6G --cpus-per-task 10 --pty bash
conda activate perly_env
# Broken pipe....better to run it using a submission script
```

```r
# Load libraries
setwd("/home/agomez/scratch/fusarium_venenatum_WD/analysis/dynGENIE3")
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
TS1 <- read.expr.matrix("new/C1R1.txt",form="rows.are.genes")
TS2 <- read.expr.matrix("new/C1R2.txt",form="rows.are.genes")
TS3 <- read.expr.matrix("new/C1R3.txt",form="rows.are.genes")
TS4 <- read.expr.matrix("new/C1R4.txt",form="rows.are.genes")
TS5 <- read.expr.matrix("new/C1R5.txt",form="rows.are.genes")

# Time
time.points <- list(TS1[1,], TS2[1,], TS3[1,], TS4[1,], TS5[1,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),], TS4[2:nrow(TS4),], TS5[2:nrow(TS5),])

# Add regulators
TFexp<-read.table("new/WT_minion_TF_expression_only.txt",header=TRUE,sep="\t")
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
ntrees <- 500
# Run the method with these settings
resall <- dynGENIE3(TS.data,time.points, regulators=TF, tree.method=tree.method, K=K, ntrees=ntrees,alpha=decay)


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

```r
# Load libraries
library(reshape2)
library(doRNG)
library(doParallel)
library(deplyr)
source ("dynGENIE3.R")

#Test works!!!

# vst blind=false expression data

rawF<-read.text("Carbon_data_vst_F.csv",header=T,sep="\t")

# NA errors!!!

# vst blind=true expression data

rawF<-read.csv("Carbon_data_vst_F.csv",header=T,sep="\t")

# Separate replicas as independent time-series experiments
C1R1<-select(rawF,1,2,7,12,17,22,27,32,37)
C1R2<-select(rawF,1,3,8,13,18,23,28,33,38)
C1R3<-select(rawF,1,4,9,14,19,24,29,34,39)
C1R4<-select(rawF,1,5,10,15,20,25,30,35,40)
C1R5<-select(rawF,1,6,11,16,21,26,31,36,41)

dat11 <-as.matrix(C1R1[-1])
dat11 <-t(dat11)
dat11<- scale(dat11, center = TRUE, scale = TRUE)
dat11<-t(dat11)
row.names(dat11)<-C1R1$X

colnames(dat11)[1]<-"0"
colnames(dat11)[2]<-"1"
colnames(dat11)[3]<-"2"
colnames(dat11)[4]<-"3"
colnames(dat11)[5]<-"4"
colnames(dat11)[6]<-"5"
colnames(dat11)[7]<-"6"
colnames(dat11)[8]<-"7"

dat12 <-as.matrix(C1R2[-1])
dat12 <-t(dat12)
dat12<- scale(dat12, center = TRUE, scale = TRUE)
dat12<-t(dat12)
row.names(dat11)<-C1R1$X

colnames(dat12)[1]<-"0"
colnames(dat12)[2]<-"1"
colnames(dat12)[3]<-"2"
colnames(dat12)[4]<-"3"
colnames(dat12)[5]<-"4"
colnames(dat12)[6]<-"5"
colnames(dat12)[7]<-"6"
colnames(dat12)[8]<-"7"



write.table(dat11, "C1R1Z.txt", sep="\t",quote = FALSE)
write.table(dat12, "C1R2Z.txt", sep="\t",quote = FALSE)


write.table(dat11, "C1R3Z.txt", sep="\t",quote = FALSE)
write.table(dat11, "C1R4Z.txt", sep="\t",quote = FALSE)
write.table(dat11, "C1R5Z.txt", sep="\t",quote = FALSE)
write.table(dat11, "C1R6Z.txt", sep="\t",quote = FALSE)
write.table(dat11, "C1R7Z.txt", sep="\t",quote = FALSE)







write.table(C1R1, "C1R1.txt", sep="\t",quote = FALSE)
write.table(C1R2, "C1R2.txt", sep="\t",quote = FALSE)
write.table(C1R3, "C1R3.txt", sep="\t",quote = FALSE)
write.table(C1R4, "C1R4.txt", sep="\t",quote = FALSE)
write.table(C1R5, "C1R5.txt", sep="\t",quote = FALSE)

# Edit headers before load matrices into R with read.expr.matrix function (from dynGenie3)
TS1 <- read.expr.matrix("C1R1.txt",form="rows.are.genes")
TS2 <- read.expr.matrix("C1R2.txt",form="rows.are.genes")
TS3 <- read.expr.matrix("C1R3.txt",form="rows.are.genes")
TS4 <- read.expr.matrix("fusarium_venenatum_WD/analysis/dynGENIE3/CarbonRNAseq/v3.0/C1R4.txt",form="rows.are.genes",stringsAsFactors = TRUE)
TS5 <- read.expr.matrix("fusarium_venenatum_WD/analysis/dynGENIE3/CarbonRNAseq/v3.0/C1R5.txt",form="rows.are.genes",stringsAsFactors = TRUE)

# Time
time.points <- list(TS1[1,], TS2[1,], TS3[1,], TS4[1,], TS5[1,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),], TS4[2:nrow(TS4),], TS5[2:nrow(TS5),])

# Add regulators
TFexp<-read.table("fusarium_venenatum_WD/analysis/dynGENIE3/CarbonRNAseq/v3.0/WT_minion_TF_expression_only.txt",header=TRUE,sep="\t")
TF<-TFexp[,1]

set.seed(123)
# Use the Extra-Trees as tree-based method
tree.method <- "ET"
# Number of randomly chosen candidate regulators at each node of a tree
K <- "all"
# Number of trees per ensemble
ntrees <- 500
# Run the method with these settings
resall <- dynGENIE3(TS.data,time.points, regulators=TF, tree.method=tree.method, K=K, ntrees=ntrees)







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

```



sbatch /mnt/shared/scratch/agomez/apps/git_repos/bioinformatics_tools/Co-expression_analysis/dynGENIE3.sh CarbonRNAseq/v4.1/C1R1.txt CarbonRNAseq/v4.1/C1R2.txt CarbonRNAseq/v4.1/C1R3.txt CarbonRNAseq/v4.1/C1R4.txt CarbonRNAseq/v4.1/C1R5.txt CarbonRNAseq/v4.1/WT_minion_TF_expression_only.txt CarbonRNAseq/v4.1



## Aracne-AP

```
Notes for Aracne.
Aracne needs a large number of samples to infer networks. Therefore, use the whole expression matrix as input (don't collapse pseudoreps!)
Create one network per condition to compare with other methods. 
I recommend vst transformed data using blind=FALSE to account for differences in expression by terms in the design (e.g. media, cultivar....)
```

Condition 1 - Glucose High

```bash
# Calculate a threshold for Mutual Information
java -Xmx5G -jar  ../../../../apps/prog/ARACNe-AP/dist/aracne.jar -e C1_all.txt -o C1_results -t WT_minion_TF_expression_only_head.txt --pvalue 1E-8 --seed 1 --calculateThreshold
# Run ARACNe on bootstraps of the input matrix
for i in {1..100}
do
java -Xmx5G -jar  ../../../../apps/prog/ARACNe-AP/dist/aracne.jar -e C1_all.txt -o C1_results -t WT_minion_TF_expression_only_head.txt --pvalue 1E-8 --seed $i
done
# Consolidate, i.e. combine the bootstraps into a final network file
java -Xmx5G -jar  ../../../../apps/prog/ARACNe-AP/dist/aracne.jar -o C1_results --consolidate
```

Condition 2 - Sucrose High

```bash
# Calculate a threshold for Mutual Information
java -Xmx5G -jar  ../../../../apps/prog/ARACNe-AP/dist/aracne.jar -e C2_all.txt -o C2_results -t WT_minion_TF_expression_only_head.txt --pvalue 1E-8 --seed 1 --calculateThreshold
# Run ARACNe on bootstraps of the input matrix
for i in {1..100}
do
java -Xmx5G -jar  ../../../../apps/prog/ARACNe-AP/dist/aracne.jar -e C2_all.txt -o C2_results -t WT_minion_TF_expression_only_head.txt --pvalue 1E-8 --seed $i
done
# Consolidate, i.e. combine the bootstraps into a final network file
java -Xmx5G -jar  ../../../../apps/prog/ARACNe-AP/dist/aracne.jar -o C2_results --consolidate
```

Condition 3 - Glucose Low

```bash
# Calculate a threshold for Mutual Information
java -Xmx5G -jar  ../../../../apps/prog/ARACNe-AP/dist/aracne.jar -e C3_all.txt -o C3_results -t WT_minion_TF_expression_only_head.txt --pvalue 1E-8 --seed 1 --calculateThreshold
# Run ARACNe on bootstraps of the input matrix
for i in {1..100}
do
java -Xmx5G -jar  ../../../../apps/prog/ARACNe-AP/dist/aracne.jar -e C3_all.txt -o C3_results -t WT_minion_TF_expression_only_head.txt --pvalue 1E-8 --seed $i
done
# Consolidate, i.e. combine the bootstraps into a final network file
java -Xmx5G -jar  ../../../../apps/prog/ARACNe-AP/dist/aracne.jar -o C3_results --consolidate
```

Condition 4 - Sucrose Low

```bash
# Calculate a threshold for Mutual Information
java -Xmx5G -jar  ../../../../apps/prog/ARACNe-AP/dist/aracne.jar -e C4_all.txt -o C4_results -t WT_minion_TF_expression_only_head.txt --pvalue 1E-8 --seed 1 --calculateThreshold
# Run ARACNe on bootstraps of the input matrix
for i in {1..100}
do
java -Xmx5G -jar  ../../../../apps/prog/ARACNe-AP/dist/aracne.jar -e C4_all.txt -o C4_results -t WT_minion_TF_expression_only_head.txt --pvalue 1E-8 --seed $i
done
# Consolidate, i.e. combine the bootstraps into a final network file
java -Xmx5G -jar  ../../../../apps/prog/ARACNe-AP/dist/aracne.jar -o C4_results --consolidate
```