#!/usr/bin/Rscript

library(optparse)
library(Biostrings)
library(SADEG, lib.loc = "/home/bonthas/R/libs/")
source("/home/armita/git_repos/emr_repos/scripts/fusarium_venenatum/codon_usage/functions2.R")

opt_list <- list(
    make_option("--cds", type="character", help="fasta file containing cds data"),
    make_option("--low", type="character", help="IDs of lowly expressed genes"),
    make_option("--high", type="character", help="IDs of highly expressed genes"),
  )
opt = parse_args(OptionParser(option_list=opt_list))


genes <- readDNAStringSet(opt$cds)
df_low <- read.delim(file=opt$low, header=F, sep="\t")
df_high <- read.delim(file=opt$high, header=F, sep="\t")

cds <- '/home/groups/harrisonlab/project_files/fusarium_venenatum/gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.cds.fasta'
genes <- readDNAStringSet(cds)
low <- 'analysis/codon_usage/fpkm_1-10/quorn_fpkm_1-10_IDs_filtered.txt'
df_low <- read.delim(file=low, header=T, sep="\t")
high <- 'analysis/codon_usage/fpkm_148/quorn_fpkm_148_IDs_filtered.txt'
df_high <- read.delim(file=high, header=T, sep="\t")

# Load libraries
#
# # Calculate average FPKM per condition
# RH <- cbind(df_fpkm[,c(1,2)],
#       rowMeans(df_fpkm[,c(3:5)]),
#       rowMeans(df_fpkm[,c(6:8)]),
#       rowMeans(df_fpkm[,c(9:11)]),
#       rowMeans(df_fpkm[,c(12:14)]),
#       rowMeans(df_fpkm[,c(15:17)]),
#       rowMeans(df_fpkm[,c(18:20)]),
#       rowMeans(df_fpkm[,c(21:23)]),
#       rowMeans(df_fpkm[,c(24:26)])
# )
# colnames(RH) <- c("GeneID","Length","RH1_mean","RH2_mean","RH3_mean","RH4_mean","RH5_mean","RH6_mean","RH7_mean","RH8_mean")

# Extract IDs from FPKM table
# and splits them into two lists of
# high and low expressed genes
# ids <- extractID(RH)

# substrRight <- function(x, n){
#   substr(x, 1, nchar(x)-n)
# }

# names(genes) <- substrRight(names(genes), 3)

ids <- list()
ids[[1]] <- df_low[,1]
ids[[2]] <- df_high[,1]

# Extract CDS sequences
# Seq <- extractSeq(ids,genes)
low_seq <- genes[ids[[1]]]
high_seq <- genes[ids[[2]]]
Seq <- list()
Seq[1] <- low_seq
Seq[2] <- high_seq

# Calculate RSCU values


matL <- lapply(Seq[[1]], function(x) {SADEG.RSCU(x)})
matH <- lapply(Seq[[2]], function(x) {SADEG.RSCU(x)})

# Process the extracted RSCU values and calculate deltaRSCU for both highly and weakly expressing genes.

# codonw[,1] <- gsub("U", "T", codonw[,1])



matH1 <- matrix(round(unlist(matH),3),ncol=64,byrow=T)
colnames(matH1) <- names(matH[[1]])
rownames(matH1) <- names(matH)
# matH1 <- matH1[,-c(38,55,62,63,64)]
matL1 <- matrix(round(unlist(matL),3),ncol=64,byrow=T)
colnames(matL1) <- names(matL[[1]])
rownames(matL1) <- names(matL)
# matL1 <- matL1[,-c(38,55,62,63,64)]


# deltaRSCU(testH=matH1,testL=matL1)
#
# deltaRSCU <- function (testH,testL,codonw) {
  var <- list()
	tmp <- list()
	# for each codon
  testL=matL1
  testH=matH1
	for(i in 1:dim(testL)[2]) {
		tmpL <- as.numeric(testL[1:dim(testL)[1],i])
		tmpH <- as.numeric(testH[1:dim(testH)[1],i])
		# tmp[[i]] <- pairwise.t.test(tmpH,tmpL,p.adjust="bonferroni")
    # Test to see if variance equal between treatments
    var[[i]] <- var.test(tmpH, tmpL)
		tmp[[i]] <- t.test(tmpH,tmpL,var.equal = F)
	}
	print(tmp[[1]])

	#make matrix 4 wide and 59 long put mean of x mean of y and p.value along with the rowname as a codon
	# Also show
	# fusTest <- matrix(data = NA, nrow = 59, ncol = 5, byrow = FALSE,dimnames = NULL)
  fusTest <- matrix(data = NA, nrow = 64, ncol = 5, byrow = FALSE,dimnames = NULL)


	# for (i in 1:59) {
	for (i in 1:64) {
		fusTest[i,1]= as.numeric(tmp[[i]]$estimate[1])
		fusTest[i,2]=as.numeric(tmp[[i]]$estimate[2])
		fusTest[i,3]=as.numeric(tmp[[i]]$p.value)
    fusTest[i,4]=as.numeric(tmp[[i]]$p.value * 64)
	}

	rownames(fusTest) <- colnames(testH)
	colnames(fusTest) <- c('high mean', 'low mean', 'p-val', 'p-adj', 'category')
	print(head(fusTest))

	# rownames(fusTest)=codon

	# print(codon)

	print(dim(fusTest))

  # fusTest[,1] <- as.numeric(fusTest[,1])
  # fusTest[,2] <- as.numeric(fusTest[,2])
  # fusTest[,3] <- as.numeric(fusTest[,3])
  # fusTest[,4] <- as.(fusTest[,4])

	# categroise into reppresented in high, low or non-sig
	# for (i in 1:59) {
  for (i in 1:64) {
		if (fusTest[i,1]>fusTest[i,2] & fusTest[i,4]<1e-2){
			fusTest[i,5]=3
		}
		if (fusTest[i,1]<fusTest[i,2] & fusTest[i,4]<1e-2){
			fusTest[i,5]=1
		}
		if (fusTest[i,1]<fusTest[i,2] & fusTest[i,4]>1e-2){
			fusTest[i,5]=2
		}
		if (fusTest[i,1]>fusTest[i,2] & fusTest[i,4]>1e-2){
			fusTest[i,5]=2
		}
	}
	print(head(fusTest))
  write.table(data.frame(fusTest), file = "fus_test.txt",append = FALSE, sep = "\t",row.names = TRUE, col.names = TRUE, quote=FALSE)

# 	arrangedfus = matrix(data = NA, nrow = 64, ncol = 2, byrow = FALSE, dimnames = NULL)
#
# 	for (i in 1:64) {
# 		if ( length(which (rownames(fusTest)==codonw[i,1]) )<1  ) {
# 			arrangedfus[i,1]=codonw[i,]
# 			arrangedfus[i,2]=2
# 		}
# 		else {
# 			arrangedfus[i,1]=rownames(fusTest)[which (rownames(fusTest)==codonw[i,])]
# 			arrangedfus[i,2]=as.numeric(fusTest[ which (rownames(fusTest)==codonw[i,]),4])
# 		}
# 	}
# 	print(arrangedfus)
# 	fusCOA =matrix(as.numeric(arrangedfus[,2]), nrow = 4, ncol = 16, byrow = TRUE, dimnames = NULL)
# 	print("Writing results to files!")
#
# 	#write.table(data.frame(arrangedfus), file = "fus_arr_new.txt",append = FALSE, sep = ",",row.names = FALSE, col.names = FALSE,qmethod = "escape")
# 	#write.table(fusCOA, file = "fus_new.coa",append = FALSE, sep = ",",row.names = FALSE, col.names = FALSE)
#
# 	return(fusCOA)
# }


# processMat(matH=matH,matL=matL, codonw=codonw)




# For codon bias per media
# processMat2(matH=matH,matL=matL,codonw=codonw)

# End of program
