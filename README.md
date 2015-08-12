# fusarium_venenatum
Bioinformatic analysis of fusarium venenatum genomes

All work was carried out in the directory:

Commands used during analysis of the neonectria_galligena genome. Note - all this work was performed in the directory:
```bash
cd /home/groups/harrisonlab/raw_data/raw_seq/fusarium/fusarium_venenatum
```

The following is a summary of the work presented in this Readme:
Data organisation:
  * Preparing data  
Draft Genome assembly
  * Data qc
  * Genome assembly
  * Repeatmasking
  * Gene prediction
  * Functional annotation
Genome analysis
  * Homology between predicted genes & published effectors


#Data organisation

Data was copied from the raw_data repository to a local directory for assembly
and annotation.

```

```


#Data qc

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:

```bash

```

Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf


```bash
  qsub /home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc/rna_qc_fastq-mcf.sh raw_dna/paired/F.venenatum/strain1/F/Fv_F_appended.fastq.gz raw_dna/paired/F.venenatum/strain1/R/Fv_R_appended.fastq.gz /home/armita/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa dna
  mv qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq.gz
  mv qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq.gz
  gunzip qc_dna/paired/F.venenatum/strain1/*/*.gz
  count_nucl.pl -i qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq -i qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq -g 60
```
  The estimated genome size is: 60000000 bp
  The input file is: qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq
  Results for: qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq
  Within this file of 2574751927 bp there were 11997418 fastq sequences
  of these 0 lines were empty.
  The input file is: qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq
  Results for: qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq
  Within this file of 2517304042 bp there were 11997418 fastq sequences
  of these 0 lines were empty.
  Total results:
  There are a total of 5092055969 nucleotides in this file.
  This equates to an estimated genome coverage of 84.87 .

Data quality was visualised once again following trimming:

```bash

```


kmer counting was performed using kmc.
This allowed estimation of sequencing depth and total genome size:

```bash
  qsub /home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc/kmer_counting.sh qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq qc_dna/paired/F.venenatum/strain1/kmer_count
```

** Estimated Genome Size is: **

** Esimated Coverage is: **

#Assembly
Assembly was performed using: Velvet / Abyss / Spades
```bash
  qsub /home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/velvet/submit_velvet_range.sh 35 65 2 qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq 60 exp_cov min_cov 600
  gzip qc_dna/paired/F.venenatum/strain1/*/*.gz
  mv qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq.gz
  mv qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq.gz
```

```bash
  F_Read=qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq
  R_Read=qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
  Outdir=assembly/spades/F.venenatum/strain1
  qsub $ProgDir/submit_SPAdes.sh $F_Read $R_Read $Outdir correct
```

Assemblies were summarised to allow the best assembly to be determined by eye.

** Assembly stats are:
  * Assembly size:
  * N50:
  * N80:
  * N20:
  * Longest contig:
  **

The assembled contigs were filtered to remove all contigs shorter than 1kb from
the assembly. This was done using the following commands:

```

```  



# Repeatmasking

Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The best assembly was used to perform repeatmasking

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
  BestAss=<PATH_TO_BEST_ASSEMBLY.fa>
  qsub $ProgDir/rep_modeling.sh $BestAss
  qsub $ProgDir/transposonPSI.sh $BestAss
```

** % bases maked by repeatmasker: **

** % bases masked by transposon psi: **


# Gene Prediction
Gene prediction followed two steps:
Pre-gene prediction - Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
Gene models were used to predict genes in the Neonectria genome. This used results from CEGMA as hints for gene models.

## Pre-gene prediction
Quality of genome assemblies was assessed by looking for the gene space in the assemblies.

```bash

```
The cegma completeness report gave an indication of the number of genes core
eukaryotic genes were present:
** Number of cegma genes present and complete: **
** Number of cegma genes present and partial: **

##Gene prediction

Gene prediction was performed for the neonectria genome.
CEGMA genes could be used as hints for the location of CDS.

For the moment we shall just use the gene model trained to F. gramminearum.
This model is from a closely related organism that is also plant pathogen.

```bash

```

** Number of genes predicted: **
The difference in drunning conditions between ERM and Nz script were assessed by
running assembly on the longest assembled contig.

```bash

```

#Functional annotation

Interproscan was used to give gene models functional annotations. Annotation was
 run using the commands below:

Note: This is a long-running script. As such, these commands were run using
 'screen' to allow jobs to be submitted and monitored in the background.
 This allows the session to be disconnected and reconnected over time.

Screen ouput detailing the progress of submission of interporscan jobs was
redirected to a temporary output file named interproscan_submission.log .


```bash

```

Following interproscan annotation split files were combined using the following commands:

```bash

```

#Genomic analysis
The first analysis was based upon BLAST searches for genes known to be involved in toxin production


##Genes with homology to PHIbase
Predicted gene models were searched against the PHIbase database using tBLASTx.

```bash

```

Top BLAST hits were used to annotate gene models.

```bash

```

following blasting PHIbase to the genome, the hits were filtered by effect on
virulence.
