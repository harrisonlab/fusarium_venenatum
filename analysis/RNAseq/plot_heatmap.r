#!/usr/bin/Rscript

annotTab_df <- read.delim("/Volumes/GGB/Harrison/Projects/BBSRC IPA QUORN C300045/Science/Obj1 Reference generation/2019/WT_annotation_ncbi_expression.tsv")
View(annotTab_df)

#install.packages("gplots")
library("gplots")
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)
setwd("~/Downloads/Fv/heatmaps")

#-----
# Heatmap plotting function
#-----


# heatmap.2(
#   rxlr_mat,
#   trace="none",
#   col=r       # use on color palette defined earlier
# )
# # This plots the heamap but we dont have control of the order of the columns

plot_heatmap <- function(matrix, prefix) {
  nrows <- dim(matrix)[1]
  png(paste(prefix, "_heatmap.png", sep=""),    # create PNG for the heat map
   width = 5*300,        # 5 x 300 pixels
   # height = (nrows*32)+400+200,
   # height = 15*300,
   height = 10*300,
   res = 300,            # 300 pixels per inch
   pointsize = 12)       # smaller font size

   heatmap.2(
     matrix,
     trace="none",
     # Colv=FALSE, # no col dendrogram and columns ordered by input
     dendrogram="row", # no col dendrogram and columns ordered by input
  #   lmat=rbind( c(0, 3, 4), c(2,1,1 ) ), lwid=c(1.5, 4, 2 ), # Change position and size of legend
     cexRow=0.75,
     cexCol=1,
    lhei = c(1,5), # Adjust relative height of key to heatmap
    # lhei = c(1,8),
    density.info="none", # remove density info from key
      margins = c(6, 6),
     col=r       # use on color palette defined earlier
   )

  dev.off()               # close the PNG device
}
#---
# End plot_heatmap function
#---


#---
# Plot Transcription factor expression
#---

tf_df <- subset(annotTab_df, TF!="" & DEG!="", select=c(1,49:56))
tf_mat <- as.matrix(tf_df[c(3,4,5,6,7,8,9)])
rownames(tf_mat) <- tf_df[,1]
colnames(tf_mat) <- c('02793', 'F55', '10170', 'MWT', 'MLO', 'MKO', 'TJ')

plot_heatmap(tf_mat, "transcription_factor")


#---
# Plot iterate through Secmet families plotting expression
#---
plot_level_all <- function(element){
  func_df <- subset(annotTab_df, Cluster_ID==element, select=c(1,49:56))
  func_mat <- as.matrix(func_df[c(3,4,5,6,7,8,9)])
  func_df$gene_id_sig <-  paste(func_df$gene_id, gsub(".+", '\\*', func_df$DEG), sep="")
  rownames(func_mat) <- func_df$gene_id_sig
  colnames(func_mat) <- c('02793', 'F55', '10170', 'MWT', 'MLO', 'MKO', 'TJ')
  if (dim(func_mat)[1] > 2){
    element <- gsub(":", "-", element)
    plot_heatmap(func_mat, element)
  }
}

plot_level_sig <- function(element){
  func_df <- subset(annotTab_df, Cluster_ID==element & DEG!="", select=c(1,49:56))
  func_mat <- as.matrix(func_df[c(3,4,5,6,7,8,9)])
  rownames(func_mat) <- func_df[,1]
  colnames(func_mat) <- c('02793', 'F55', '10170', 'MWT', 'MLO', 'MKO', 'TJ')
  if (dim(func_mat)[1] > 2){
    element <- gsub(":", "-", element)
    plot_heatmap(func_mat, element)
  }
}

x <- levels(annotTab_df$Cluster_ID)

lapply(x[-1], plot_level_all)
# lapply(x[-1], plot_level_sig)
