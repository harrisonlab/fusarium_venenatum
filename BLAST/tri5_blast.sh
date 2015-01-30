#!/bin/bash
# Commands used to perform blast searching in the F. venenatum genome for the tri5 gene.
mkdir -p ../../../../project_files/fusarium_venenatum/raw_dna/paired/F.venenatum/strain1/F/
mkdir -p ../../../../project_files/fusarium_venenatum/raw_dna/paired/F.venenatum/strain1/R/
cp fungus1_S1_L001_R1_001.fastq.gz ../../../../project_files/fusarium_venenatum/raw_dna/paired/F.venenatum/strain1/F/.
cp fungus1_S1_L001_R2_001.fastq.gz ../../../../project_files/fusarium_venenatum/raw_dna/paired/F.venenatum/strain1/R/.
cp fungus2_S1_L001_R1_001.fastq.gz ../../../../project_files/fusarium_venenatum/raw_dna/paired/F.venenatum/strain1/F/.
cp fungus2_S1_L001_R2_001.fastq.gz ../../../../project_files/fusarium_venenatum/raw_dna/paired/F.venenatum/strain1/R/.
cp *.zip ../../../../project_files/fusarium_venenatum/assembly/other_groups/F.venenatum/strain1/.

#Download Tri5 gene from genbank

# BLAST against already assembled genome
unzip assembly/other_groups/F.venenatum/strain1/1_S1_L001_R1_001\ \(paired\)\ trimmed\ \(paired\)\ assembly.fa.zipmv 1_S1_L001_R1_001\ \(paired\)\ trimmed\ \(paired\)\ assembly.fa assembly/other_groups/F.venenatum/strain1/Fv_assembly.faqsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/genbank/g.zea_Tri5.fasta protein assembly/other_groups/F.venenatum/strain1/Fv_assembly.fa 
less analysis/blast_homology/other_groups/F.venenatum/F.venenatum_g.zea_Tri5.fasta_homologs.csv

# Perform our own genome assembly and analysis
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc/rna_qc_fastq-mcf.sh raw_dna/paired/F.venenatum/strain1/F/Fv_F_appended.fastq.gz raw_dna/paired/F.venenatum/strain1/R/Fv_R_appended.fastq.gz /home/armita/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa