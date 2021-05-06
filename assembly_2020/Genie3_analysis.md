# Genie3 analysis

### Analysis of one modele from WGCNA

```r
# Convert TFs to genes ID
less analysis/transcription_factors/F.venenatum/WT_minion/WT_minion_TF_gene_headers.txt | sed "s/\..*//" > WT_minion_TF_gene_only_headers.txt

# Modele with TRI5 gene cluster
T1<-read.table("analysis/WGCNA/merged_modules/navajowhite2_cleaned.txt",header=T,sep="\t")
# All TFs
T2<-read.table("analysis/transcription_factors/F.venenatum/WT_minion/WT_minion_TF_gene_only_headers.txt",header=T)
# TFs in the module for regulation
T3<-merge(T1,T2, by.x="ID",by.y="ID",all.x=FALSE,all.y=FALSE)
# All data
T4<-read.table("analysis/WGCNA/vst1_corrected.txt",header=T,sep="\t")
# Extract expression data of genes in the module
T5<-merge(T1,T4, by.x="ID",by.y="ID",all.x=FALSE,all.y=FALSE)
write.table(T5, "TRI5module.txt", sep="\t")
# Remove ID from T5
T6<-read.table("TRI5module.txt",header=T)

# All genes candidate regulator. Need a matrix format
weightMat <- GENIE3(as.matrix(T6))
# Check input
dim(weightMat)
weightMat[1:5,1:5]

# Define regulators
regulators <- c("g4106","g4134","g5413","g6132")
# Genie3 using TFs as regulators
weightMat <- GENIE3(as.matrix(T6), regulators=regulators)

# Check input
dim(weightMat)
[1]  4 59
weightMat[1:4,1:4]

# Transform results to network

# All genes in the module as candidates

library(GENIE3)
library(igraph)
library(RCy3)
library(Rgraphviz)

weightMat <- GENIE3(as.matrix(T6))
set.seed(weightMat)

linkList <- getLinkList(weightMat)

edge_listsi <- linkList
# Build graph from dataframe
Gsi <- graph.data.frame(edge_listsi,directed = F)
# Convert graph to adjacency matrix
Asi <- get.adjacency(Gsi,sparse = F,attr = "weight",type = "both")
# Build adjacency graph
g_arasi <- graph.adjacency(Asi,mode = "undirected",weighted = T)
# Create igraph
g.cyto <- igraph.to.graphNEL(g_arasi)

cw = createNetworkFromGraph(graph=g.cyto)
#displayGraph (cw)

#######

# TFs regulators 

regulators <- c("g4106","g4134","g5413","g6132")
weightMat <- GENIE3(as.matrix(T6), regulators=regulators)

linkList <- getLinkList(weightMat)

edge_listsi <- linkList
# Build graph from dataframe
Gsi <- graph.data.frame(edge_listsi,directed = F)
# Convert graph to adjacency matrix
Asi <- get.adjacency(Gsi,sparse = F,attr = "weight",type = "both")
# Build adjacency graph
g_arasi <- graph.adjacency(Asi,mode = "undirected",weighted = T)
# Create igraph
g.cyto <- igraph.to.graphNEL(g_arasi)

cw = createNetworkFromGraph(graph=g.cyto)
displayGraph (cw)


#####

# Added for reproducibility of results
set.seed(123)
regulators <- c("g4106","g4134","g5413","g6132")
weightMat <- GENIE3(as.matrix(T6), regulators=regulators)
linkList <- getLinkList(weightMat, threshold=0.4)

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
displayGraph (cw)

```


# Analysis with all the genes

```r

less analysis/transcription_factors/F.venenatum/WT_minion/WT_minion_TF_gene_only_headers2.txt | uniq > analysis/transcription_factors/F.venenatum/WT_minion/WT_minion_TF_gene_only_headers3.txt

# All data without ID
T8<-read.table("analysis/WGCNA/vst1_4genie3.txt",header=T)
# All TFs no duplicate no ID
T9<-read.table("analysis/transcription_factors/F.venenatum/WT_minion/WT_minion_TF_gene_only_headers3.txt",header=TRUE,sep="\t")
# All TFs no duplicate
T11<-read.table("analysis/transcription_factors/F.venenatum/WT_minion/WT_minion_TF_gene_only_headers4.txt",header=TRUE,sep="\t")

# TFs genes only with all expression data
T10<-merge(T4,T11, by.x="ID",by.y="ID",all.x=FALSE,all.y=FALSE)
# Only TFs that are in the dataset. For regulators
T15<-T10[,1]
# Write table, remove ID and input again (not needed)
write.table(T10, "analysis/transcription_factors/F.venenatum/WT_minion/WT_minion_TF_gene_only_headers5.txt", sep="\t")
T12<-read.table("analysis/transcription_factors/F.venenatum/WT_minion/WT_minion_TF_gene_only_headers5.txt",header=TRUE,sep="\t")

T15<-T12[,1]
# Regulatos final 
df <-as.character(unlist(T15))
regulators <- df
# Final call. All data vs TFs with some expression
weightMat <- GENIE3(as.matrix(T8), regulators=regulators)

linkList <- getLinkList(weightMat, threshold=0.05)

edge_listsi <- linkList
# Build graph from dataframe
Gsi <- graph.data.frame(edge_listsi,directed = F)
# Convert graph to adjacency matrix
Asi <- get.adjacency(Gsi,sparse = F,attr = "weight",type = "both")
# Build adjacency graph
g_arasi <- graph.adjacency(Asi,mode = "undirected",weighted = T)
# Create igraph
g.cyto <- igraph.to.graphNEL(g_arasi)

cw = createNetworkFromGraph(graph=g.cyto)
```