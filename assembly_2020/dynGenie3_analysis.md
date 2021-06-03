cd analysis/dynGENIE3

```bash

/projects/software/R-3.6.1/bin/R
```

```R
setwd("/projects/fusarium_venenatum/analysis/dynGENIE3")

library(GENIE3)
library(igraph)
library(RCy3)
library(Rgraphviz)

library ("reshape2")
library ("doRNG")
library ("doParallel")

source ("dynGENIE3.R")

D1<-read.table("C1_v3.txt",header=T)
D2<-read.table("C2_v3.txt",header=T)
D3<-read.table("C3_v3.txt",header=T)
D4<-read.table("C4_v3.txt",header=T)
D5<-read.table("SecmetTFs_headers.txt",header=T)

P1<-merge(D1,D5, by.x="ID",by.y="ID",all.y=TRUE)
P2<-merge(D2,D5, by.x="ID",by.y="ID",all.y=TRUE)
P3<-merge(D3,D5, by.x="ID",by.y="ID",all.y=TRUE)
P4<-merge(D4,D5, by.x="ID",by.y="ID",all.y=TRUE)

write.table(P1, "C1_v2.txt", sep="\t",quote = FALSE)
write.table(P2, "C2_v2.txt", sep="\t",quote = FALSE)
write.table(P3, "C3_v2.txt", sep="\t",quote = FALSE)
write.table(P4, "C4_v2.txt", sep="\t",quote = FALSE)

TS1 <- read.expr.matrix("C1_v2.txt",form="rows.are.samples")
TS2 <- read.expr.matrix("C2_v2.txt",form="rows.are.samples")
TS3 <- read.expr.matrix("C3_v2.txt",form="rows.are.samples")
TS4 <- read.expr.matrix("C4_v2.txt",form="rows.are.samples")

time.points <- list(TS1[1,], TS2[1,], TS3[1,], TS4[1,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),], TS4[2:nrow(TS4),])

# Add regulators
T12<-read.table("WT_minion_TF_gene_only_headers5.txt",header=TRUE,sep="\t")
T15<-T12[,1]
# Regulatos final 
df <-as.character(unlist(T15))
regulators <- df

# Run dynGENIE3
res <- dynGENIE3(TS.data,time.points, regulators=regulators)


link.list <- get.link.list(res$weight.matrix)
head(link.list)

link.list <- get.link.list(res$weight.matrix, threshold=0.01)
write.table(link.list, "dyGENIE3_001.txt", sep="\t")

link.list <- get.link.list(res$weight.matrix, threshold=0.02)


linkList2 <- getLinkList(weightMat, threshold=0.05)

edge_listsi <- link.list
# Build graph from dataframe
Gsi <- graph.data.frame(edge_listsi,directed = F)
# Convert graph to adjacency matrix
Asi <- get.adjacency(Gsi,sparse = F,attr = "weight",type = "both")
# Build adjacency graph
g_arasi <- graph.adjacency(Asi,mode = "undirected",weighted = T)
# Create igraph
g.cyto <- igraph.to.graphNEL(g_arasi)

cw = createNetworkFromGraph(graph=g.cyto)


# High nitrogen only

time.points <- list(TS1[1,], TS2[1,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),])

time.points <- list(TS2[1,])
TS.data <- list(TS2[2:nrow(TS2),])

set.seed(123)
# Use the Extra-Trees as tree-based method
tree.method <- "RF"
# Number of randomly chosen candidate regulators at each node of a tree
K <- "all"
# Number of trees per ensemble
ntrees <- 1000
# Run the method with these settings
res8 <- dynGENIE3(TS.data,time.points, regulators=regulators, tree.method=tree.method, K=K, ntrees=ntrees)

link.list <- get.link.list(res8$weight.matrix)
head(link.list)
link.list <- get.link.list(res8$weight.matrix, threshold=0.1)

write.table(link.list, "dyGENIE3_RF_01_1000trees.txt", sep="\t")


# TRI5 cluster only

D1<-read.table("C1_v3.txt",header=T)
D2<-read.table("C2_v3.txt",header=T)
D3<-read.table("C3_v3.txt",header=T)
D4<-read.table("C4_v3.txt",header=T)
D5<-read.table("TF_TRI.txt",header=T)

P1<-merge(D1,D5, by.x="ID",by.y="ID",all.y=TRUE)
P2<-merge(D2,D5, by.x="ID",by.y="ID",all.y=TRUE)
P3<-merge(D3,D5, by.x="ID",by.y="ID",all.y=TRUE)
P4<-merge(D4,D5, by.x="ID",by.y="ID",all.y=TRUE)

write.table(P1, "C1_v4.txt", sep="\t",quote = FALSE)
write.table(P2, "C2_v4.txt", sep="\t",quote = FALSE)
write.table(P3, "C3_v4.txt", sep="\t",quote = FALSE)
write.table(P4, "C4_v4.txt", sep="\t",quote = FALSE)


TS1 <- read.expr.matrix("C1_v5.txt",form="rows.are.samples")
TS2 <- read.expr.matrix("C2_v5.txt",form="rows.are.samples")
TS3 <- read.expr.matrix("C3_v5.txt",form="rows.are.samples")
TS4 <- read.expr.matrix("C4_v5.txt",form="rows.are.samples")

time.points <- list(TS1[1,], TS2[1,], TS3[1,], TS4[1,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),], TS4[2:nrow(TS4),])

# Add regulators
T12<-read.table("WT_minion_TF_gene_only_headers5.txt",header=TRUE,sep="\t")
T15<-T12[,1]
# Regulatos final 
df <-as.character(unlist(T15))
regulators <- df

set.seed(123)
# Use the Extra-Trees as tree-based method
tree.method <- "RF"
# Number of randomly chosen candidate regulators at each node of a tree
K <- "all"
# Number of trees per ensemble
ntrees <- 1000
# Run the method with these settings
restri <- dynGENIE3(TS.data,time.points, regulators=regulators, tree.method=tree.method, K=K, ntrees=ntrees)

link.list <- get.link.list(restri$weight.matrix)
head(link.list)

link.list <- get.link.list(restri$weight.matrix, threshold=0.1)
write.table(link.list, "dyGENIE3_RF_tri_01.txt", sep="\t")

set.seed(123)
# Use the Extra-Trees as tree-based method
tree.method <- "ET"
# Number of randomly chosen candidate regulators at each node of a tree
K <- "all"
# Number of trees per ensemble
ntrees <- 1000
# Run the method with these settings
restri2 <- dynGENIE3(TS.data,time.points, regulators=regulators, tree.method=tree.method, K=K, ntrees=ntrees)

link.list <- get.link.list(restri2$weight.matrix)
head(link.list)

link.list <- get.link.list(restri2$weight.matrix, threshold=0.05)
write.table(link.list, "dyGENIE3_ET_tri_005.txt", sep="\t")

link.list <- get.link.list(restri2$weight.matrix, threshold=0.1)

edge_listsi <- link.list
# Build graph from dataframe
Gsi <- graph.data.frame(edge_listsi,directed = F)
# Convert graph to adjacency matrix
Asi <- get.adjacency(Gsi,sparse = F,attr = "weight",type = "both")
# Build adjacency graph
g_arasi <- graph.adjacency(Asi,mode = "undirected",weighted = T)
# Create igraph
g.cyto <- igraph.to.graphNEL(g_arasi)

cw = createNetworkFromGraph(graph=g.cyto)


# plot 

Z1<-read.table("TFs.txt",header=T)
Z2<-read.table("vst1_corrected.txt",header=T)


Z3<-merge(Z1,Z2, by.x="ID",by.y="ID",all.x=TRUE)
reshaped2 <- melt(Z3, id=c("ID"), variable.name="Timepoint", value.name="vst")
write.csv(reshaped2, "TRI5_TFs.csv")

Tri <- read.table("TRI5_TFs.txt",header=T,sep="\t")

ggplot(Tri, aes(Timepoint,vst, group=ID, color=ID))+
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


# Toxins only

time.points <- list(TS1[1,], TS2[1,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),])

time.points <- list(TS2[1,])
TS.data <- list(TS2[2:nrow(TS2),])

set.seed(123)
# Use the Extra-Trees as tree-based method
tree.method <- "ET"
# Number of randomly chosen candidate regulators at each node of a tree
K <- "all"
# Number of trees per ensemble
ntrees <- 50
# Run the method with these settings
res9 <- dynGENIE3(TS.data,time.points, regulators=regulators, tree.method=tree.method, K=K, ntrees=ntrees)

link.list <- get.link.list(res9$weight.matrix)
head(link.list)
link.list <- get.link.list(res9$weight.matrix, threshold=0.03)

write.table(link.list, "dyGENIE3_RF_01_1000trees.txt", sep="\t")
