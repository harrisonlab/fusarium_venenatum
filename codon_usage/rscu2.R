#!/usr/bin/Rscript

library(optparse)
library(Biostrings)
library(SADEG, lib.loc = "/home/bonthas/R/libs/")
source("/home/armita/git_repos/emr_repos/scripts/fusarium_venenatum/codon_usage/functions2.R")

opt_list <- list(
    make_option("--cds", type="character", help="fasta file containing cds data"),
    make_option("--fpkm", type="character", help="tab seperated fpkm data, per condition"),
    make_option("--codons", type="character", help="codonw codon list")
  )
opt = parse_args(OptionParser(option_list=opt_list))

#
# genes <- readDNAStringSet(opt$cds)
# df_fpkm <- read.delim(file=opt$fpkm, header=F, sep="\t")
# codonw <- read.delim(file=opt$codons, header=F, sep="\n", quote="",stringsAsFactors = FALSE)
cds <- '/home/groups/harrisonlab/project_files/fusarium_venenatum/gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.cds.fasta'
genes <- readDNAStringSet(cds)
fpkm <- '/home/deakig/projects/quorn/DGE/Quorn_fpkm.txt'
df_fpkm <- read.delim(file=fpkm, header=T, sep="\t")
codons <- '/home/groups/harrisonlab/project_files/fusarium_venenatum/analysis/codon_usage/final_genes_appended_renamed.cds/codon.coa'
codonw <- read.delim(file=codons, header=T, sep=" ", quote="",stringsAsFactors = FALSE)

# Load libraries

# Calculate average FPKM per condition
RH <- cbind(df_fpkm[,c(1,2)],
      rowMeans(df_fpkm[,c(3:5)]),
      rowMeans(df_fpkm[,c(6:8)]),
      rowMeans(df_fpkm[,c(9:11)]),
      rowMeans(df_fpkm[,c(12:14)]),
      rowMeans(df_fpkm[,c(15:17)]),
      rowMeans(df_fpkm[,c(18:20)]),
      rowMeans(df_fpkm[,c(21:23)]),
      rowMeans(df_fpkm[,c(24:26)])
)
colnames(RH) <- c("GeneID","Length","RH1_mean","RH2_mean","RH3_mean","RH4_mean","RH5_mean","RH6_mean","RH7_mean","RH8_mean")

# Extract IDs from FPKM table
# and splits them into two lists of
# high and low expressed genes
ids <- extractID(RH)

substrRight <- function(x, n){
  substr(x, 1, nchar(x)-n)
}

names(genes) <- substrRight(names(genes), 3)

# Extract CDS sequences
Seq <- extractSeq(ids,genes)

# Calculate RSCU values
matH <- lapply(Seq[[1]], function(x) {lapply(x,SADEG.RSCU)})

matL <- lapply(Seq[[2]], function(x) {lapply(x,SADEG.RSCU)})

# Process the extracted RSCU values and calculate deltaRSCU for both highly and weakly expressing genes.

codonw[,1] <- gsub("U", "T", codonw[,1])
processMat(matH=matH,matL=matL, codonw=codonw)

# For codon bias per media
# processMat2(matH=matH,matL=matL,codonw=codonw)

# End of program
