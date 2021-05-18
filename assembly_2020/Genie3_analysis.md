# Genie3 analysis

```r
# Load libraries
library(GENIE3)
library(igraph)
library(RCy3)
library(Rgraphviz)
```

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


# Analysis with secmet genes and TFs

```r

# Modele with TRI5 gene cluster
D1<-read.table("analysis/secondary_metabolites/antismash/F.venenatum/WT_minion_VP/WT_antismash_results_secmet_genes_with_header.txt",header=T)
# All TFs
D2<-read.table("analysis/secondary_metabolites/smurf/F.venenatum/WT_minion/WT_minion_smurf_secmet_genes_withhead.txt",header=T)
# TFs in the module for regulation
D3<-merge(D1,D2, by.x="ID",by.y="ID",all.x=TRUE,all.y=TRUE)
write.table(D3, "analysis/Genie3/Secmet/WT_minion_all_secmet.txt", sep="\t",quote = FALSE)

less analysis/Genie3/Secmet/WT_minion_all_secmet.txt | cut -f1 | sed "s/\..*//" | uniq > analysis/Genie3/Secmet/combined_secmet_gene_names.txt

D4<-read.table("analysis/Genie3/Secmet/combined_secmet_gene_names.txt",header=T)
# All TFs no duplicate
D5<-read.table("analysis/transcription_factors/F.venenatum/WT_minion/WT_minion_TF_gene_only_headers4.txt",header=T)
# Secmet and TFs
D6<-merge(D4,D5,by.x="ID",by.y="ID",all.x=TRUE,all.y=TRUE)

# All data
T4<-read.table("analysis/WGCNA/vst1_corrected.txt",header=T,sep="\t")

# Extract expression data of Secmet and TFs
D7<-merge(D6,T4, by.x="ID",by.y="ID",all.x=FALSE,all.y=FALSE)
write.table(D7, "analysis/Genie3/Secmet/SecmetTFs.txt", sep="\t")

# TFs genes only with all expression data
D8<-merge(D7,D5, by.x="ID",by.y="ID",all.x=FALSE,all.y=FALSE)
D9<-D8[,1]

# Load corrected table
D10<-read.table("analysis/Genie3/Secmet/SecmetTFs.txt",header=TRUE)

# GENIE3

# Added for reproducibility of results
set.seed(123)
weightMat <- GENIE3(as.matrix(D10), regulators=D9)

#linkList <- getLinkList(weightMat)

linkList <- getLinkList(weightMat, threshold=0.03)

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

write.table(linkList, "analysis/Genie3/Secmet/Results_table.txt", sep="\t",quote = FALSE)
```

