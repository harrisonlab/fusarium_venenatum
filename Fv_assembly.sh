# Fusarium Venenatum

cd /home/groups/harrisonlab/raw_data/raw_seq/fusarium/fusarium_venenatum/

# Perform our own genome assembly and analysis
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc/rna_qc_fastq-mcf.sh raw_dna/paired/F.venenatum/strain1/F/Fv_F_appended.fastq.gz raw_dna/paired/F.venenatum/strain1/R/Fv_R_appended.fastq.gz /home/armita/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa dna
mv qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq.gz
mv qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq.gz
gunzip qc_dna/paired/F.venenatum/strain1/*/*.gz
count_nucl.pl -i qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq -i qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq -g 60
The estimated genome size is: 60000000 bp
# The input file is: qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq
# Results for: qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq
#  Within this file of 2574751927 bp there were 11997418 fastq sequences
#  of these 0 lines were empty.
# The input file is: qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq
# Results for: qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq
#  Within this file of 2517304042 bp there were 11997418 fastq sequences
#  of these 0 lines were empty.
# Total results:
#  There are a total of 5092055969 nucleotides in this file.
#  This equates to an estimated genome coverage of 84.87 .
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/velvet/submit_velvet_range.sh 35 65 2 qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq 60 exp_cov min_cov 600
gzip qc_dna/paired/F.venenatum/strain1/*/*.gz
mv qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq.gz
mv qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq.gz

qsub /home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc/kmer_counting.sh qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq qc_dna/paired/F.venenatum/strain1/kmer_count
