cd analysis/dynGENIE3

```bash

/projects/software/R-3.6.1/bin/R
```

```R
setwd("/projects/fusarium_venenatum/analysis/dynGENIE3")

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