# Dynamic regulatory networks

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

```r
# Run on Crop diversity

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



# Aracne

java -Xmx5G -jar apps/prog/ARACNe-AP/dist/aracne.jar -e fusarium_venenatum_WD/analysis/ARACNe/vtest/C1_av.txt -o fusarium_venenatum_WD/analysis/ARACNe/vtest/outputFolder --tfs fusarium_venenatum_WD/analysis/ARACNe/vtest/WT_minion_TF_expression_only_head.txt --pvalue 1E-8 --seed 1 --calculateThreshold

for i in {1..100}
do
java -Xmx5G -jar apps/prog/ARACNe-AP/dist/aracne.jar -e fusarium_venenatum_WD/analysis/ARACNe/vtest/C1_av.txt -o fusarium_venenatum_WD/analysis/ARACNe/vtest/outputFolder --tfs fusarium_venenatum_WD/analysis/ARACNe/vtest/WT_minion_TF_expression_only_head.txt --pvalue 1E-8 --seed $i
done

java -Xmx5G -jar apps/prog/ARACNe-AP/dist/aracne.jar -o fusarium_venenatum_WD/analysis/ARACNe/vtest/outputFolder --consolidate









C1R1<-select(rawF,1,2,7,12,17,22,27,32,37)
C1R2<-select(rawF,1,3,8,13,18,23,28,33,38)
C1R3<-select(rawF,1,4,9,14,19,24,29,34,39)
C1R4<-select(rawF,1,5,10,15,20,25,30,35,40)
C1R5<-select(rawF,1,6,11,16,21,26,31,36,41)

write.table(C1R1, "C1R1.txt", sep="\t",quote = FALSE)
write.table(C1R2, "C1R2.txt", sep="\t",quote = FALSE)
write.table(C1R3, "C1R3.txt", sep="\t",quote = FALSE)
write.table(C1R4, "C1R4.txt", sep="\t",quote = FALSE)
write.table(C1R5, "C1R5.txt", sep="\t",quote = FALSE)

TS1 <- read.expr.matrix("fusarium_venenatum_WD/analysis/dynGENIE3/CarbonRNAseq/v3.0/C1R1.txt",form="rows.are.genes")
TS2 <- read.expr.matrix("fusarium_venenatum_WD/analysis/dynGENIE3/CarbonRNAseq/v3.0/C1R2.txt",form="rows.are.genes")
TS3 <- read.expr.matrix("fusarium_venenatum_WD/analysis/dynGENIE3/CarbonRNAseq/v3.0/C1R3.txt",form="rows.are.genes")
TS4 <- read.expr.matrix("fusarium_venenatum_WD/analysis/dynGENIE3/CarbonRNAseq/v3.0/C1R4.txt",form="rows.are.genes")
TS5 <- read.expr.matrix("fusarium_venenatum_WD/analysis/dynGENIE3/CarbonRNAseq/v3.0/C1R5.txt",form="rows.are.genes")

time.points <- list(TS1[1,], TS2[1,], TS3[1,], TS4[1,], TS5[1,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),], TS4[2:nrow(TS4),], TS5[2:nrow(TS5),])

set.seed(123)
# Use the Extra-Trees as tree-based method
tree.method <- "ET"
# Number of randomly chosen candidate regulators at each node of a tree
K <- "all"
# Number of trees per ensemble
ntrees <- 500
# Run the method with these settings
resall <- dynGENIE3(TS.data,time.points, regulators=TF, tree.method=tree.method, K=K, ntrees=ntrees)




library(reshape2)
library(doRNG)
library(doParallel)
source ("dynGENIE3_R_C_wrapper/dynGENIE3.R")


# Edit headers before load matrices into R with read.expr.matrix function (from dynGenie3)
TS1 <- read.expr.matrix("C1R1_D.txt",form="rows.are.genes")
TS2 <- read.expr.matrix("C1R2_D.txt",form="rows.are.genes")
TS3 <- read.expr.matrix("C1R3.txt",form="rows.are.genes")
TS4 <- read.expr.matrix("C1R4.txt",form="rows.are.genes")
TS5 <- read.expr.matrix("C1R5.txt",form="rows.are.genes")



colnames(C1R1)[1]<-"gene"
colnames(C1R1)[2]<-"0"
colnames(C1R1)[3]<-"1"
colnames(C1R1)[4]<-"2"
colnames(C1R1)[5]<-"3"
colnames(C1R1)[6]<-"4"
colnames(C1R1)[7]<-"5"
colnames(C1R1)[8]<-"6"
colnames(C1R1)[9]<-"7"

colnames(C1R2)[1]<-"gene"
colnames(C1R2)[2]<-"0"
colnames(C1R2)[3]<-"1"
colnames(C1R2)[4]<-"2"
colnames(C1R2)[5]<-"3"
colnames(C1R2)[6]<-"4"
colnames(C1R2)[7]<-"5"
colnames(C1R2)[8]<-"6"
colnames(C1R2)[9]<-"7"

write.table(C1R1, "C1R1.txt", sep="\t",quote = FALSE)
write.table(C1R2, "C1R2.txt", sep="\t",quote = FALSE)
write.table(C1R3, "C1R3.txt", sep="\t",quote = FALSE)
write.table(C1R4, "C1R4.txt", sep="\t",quote = FALSE)
write.table(C1R5, "C1R5.txt", sep="\t",quote = FALSE)

#T Time
time.points <- list(TS1[1,], TS2[1,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),])

res <- dynGENIE3(TS.data,time.points, regulators=TF)


TS1 <- read.expr.matrix("test1.txt",form="rows.are.genes")
TS2 <- read.expr.matrix("test2.txt",form="rows.are.genes")