# Network diagram for common regulatory motifs between SecMet clusters.
# Commands as run from my local computer.

​
install.packages('igraph')
install.packages('network')
install.packages('sna')
install.packages('ndtv')
install.packages('visNetwork')
​
library(igraph)
library(network)
library(sna)
library(ndtv)
library(visNetwork)
library(plyr)
​
setwd("/Users/armita/Downloads/tomtom")
df2 <- read.delim("~/Downloads/tomtom/tomtom_hits_high_score.tsv", header=FALSE)
library(dplyr)
df3 <- summarise(group_by(df2,V1,V4),count =n())
library('igraph')
net <- graph_from_data_frame(d=df3, directed=T)
net
as_edgelist(net, names=T)
as_adjacency_matrix(net, attr="count")
# as_data_frame(net, what="edges")
# as_data_frame(net, what="vertices")
plot(net)

df3$V1 <- gsub('SecMet_cluster_', '', df3$V1)
df3$V4 <- gsub('SecMet_cluster_', '', df3$V4)

net <- graph_from_data_frame(d=df3, directed=F)
plot(net)
E(net)$width <- E(net)$count
plot(net,vertex.color="white", edge.color="dark red")
plot(net, vertex.label.dist=1,layout=layout_in_circle)

net <- graph_from_data_frame(d=df3, directed=T)
E(net)$width <- sqrt(E(net)$count)
#E(net)$edge.arrow.size <- sqrt(E(net)$count)

plot(net,edge.arrow.size=.4, edge.curved=.5, edge.label=df3$count, edge.label.font = 3)

pdf("Fv_SecMet_network.pdf",width=10,height=10,paper='special')
plot(net,edge.arrow.size=.4, edge.curved=.5, edge.label=df3$count, edge.label.font = 3)
dev.off()
