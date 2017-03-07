# Colletotrichum gloeosporioides
==========

Scripts used for the analysis of Colletotrichum gloeosporioides genomes
Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/colletotrichum_gloeosporioides

The following is a summary of the work presented in this Readme.

The following processes were applied to Fusarium genomes prior to analysis:
Data qc
Genome assembly
Repeatmasking
Gene prediction
Functional annotation

<!--
Analyses performed on these genomes involved BLAST searching for:

 ls contigs were identified using:
Alignment of raw reads to assembled genomes
Assembly of remaining reads
-->


#Building of directory structure


```bash
  # Oxford nanopore 07/03/17
  RawDatDir=/home/miseq_data/minion/2017/Fvenenatum/downloaded/pass
  Species=C.venenatum
  Strain=WT
  Date=07-03-17
  mkdir -p raw_dna/minion/$Species/$Strain/$Date
  for Fast5Dir in $(ls -d $RawDatDir/*); do
    poretools fastq $Fast5Dir | gzip -cf
  done > raw_dna/minion/$Species/$Strain/"$Strain"_"$Date"_fail.fastq.gz
  mv raw_dna/minion/$Species/$Strain/"$Strain"_"$Date"_fail.fastq.gz raw_dna/minion/$Species/$Strain/"$Strain"_"$Date"_pass.fastq.gz
  # poretools stats $RawDatDir/ > raw_dna/minion/$Species/$Strain/"$Strain"_"$Date"_fail.stats.txt
  # poretools hist $RawDatDir/ > raw_dna/minion/$Species/$Strain/"$Strain"_"$Date"_fail.hist
  # cat raw_dna/minion/$Species/$Strain/"$Strain"_"$Date".fastq.gz raw_dna/minion/$Species/$Strain/"$Strain"_"$Date"_fail.fastq.gz > raw_dna/minion/$Species/$Strain/"$Strain"_"$Date"_pass-fail.fastq.gz
```
