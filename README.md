# fusarium_venenatum
Bioinformatic analysis of fusarium venenatum genomes

All work was carried out in the directory:

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
Assembly was performed using: Spades
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

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  Assembly=assembly/spades/F.venenatum/strain1/filtered_contigs/contigs_min_500bp.fasta
  OutDir=assembly/spades/F.venenatum/strain1/filtered_contigs
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
```

As SPADes was run with the option to autodetect a minimum coverage the assembly
was assessed to identify the coverage of assembled contigs. This was done using
the following command:

```
  BestAss=assembly/spades/F.venenatum/strain1/filtered_contigs/contigs_min_500bp.fasta
  cat $BestAss | grep '>' | cut -f6 -d'_' | sort -n | cut -f1 -d '.' | sort -n | uniq -c | less
```

This returned the following results:
```
  633 2
   85 3
    4 4
    2 5
    2 6
    2 7
    1 8
    1 11
    2 13
    1 14
    3 15
    2 16
    2 18
    1 19
    1 20
    2 23
    1 27
    1 28
    2 33
    3 36
    2 38
    1 39
    4 40
    3 41
    7 42
    6 43
   21 44
   46 45
   50 46
   14 47
    1 48
    2 49
    1 53
    1 54
    1 55
    1 57
    1 58
    1 60
    1 63
    1 68
    1 80
    1 87
    1 92
    1 104
    1 105
    1 108
    1 109
    1 137
    1 170
    1 354
    1 503
    1 506
    1 566
    1 791
    1 903
    1 987
    1 1025
    1 1028
    1 1119
    1 1136
    1 1179
    1 1263
    1 1484
    1 1553
    1 4306
    1 5188
```
From this it was determined that SPades could not be trusted to set its own
minimum threshold for coverage.

In future an option will be be used to set a coverage for spades.

in the meantime contigs with a coverage lowere than 10 were filtered out using
the following commands:

```
  Headers=assembly/spades/F.venenatum/strain1/filtered_contigs/contigs_min_500bp_10x_headers.txt
  FastaMinCov=assembly/spades/F.venenatum/strain1/filtered_contigs/contigs_min_500bp_10x.fasta
  cat $BestAss | grep '>' | grep -E -v 'cov_.\..*_' > $Headers

  cat $BestAss | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep -A1 -f $Headers | grep -v -E '^\-\-' > $FastaMinCov
```

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  Assembly=assembly/spades/F.venenatum/strain1/filtered_contigs/contigs_min_500bp_10x.fasta
  OutDir=assembly/spades/F.venenatum/strain1/filtered_contigs/contigs_min_500bp_10x
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
```

Results were as follows:
```
  All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

  Assembly                   contigs_min_500bp_10x  contigs_min_500bp_10x broken
  # contigs (>= 0 bp)        213                    213
  # contigs (>= 1000 bp)     196                    196
  Total length (>= 0 bp)     37650527               37650527
  Total length (>= 1000 bp)  37638009               37638009
  # contigs                  212                    212
  Largest contig             1154964                1154964
  Total length               37650028               37650028
  GC (%)                     47.62                  47.62
  N50                        392022                 392022
  N75                        233930                 233930
  L50                        32                     32
  L75                        62                     62
  # N's per 100 kbp          0.00                   0.00
```

Contigs were renamed in accordance with ncbi recomendations.

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  for Assembly in $(ls assembly/spades/F.venenatum/strain1/filtered_contigs/contigs_min_500bp_10x.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

# Repeatmasking

Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The best assembly was used to perform repeatmasking

```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/repeat_masking
  BestAss=assembly/spades/F.venenatum/strain1/filtered_contigs/contigs_min_500bp_renamed.fasta
  qsub $ProgDir/rep_modeling.sh $BestAss
  qsub $ProgDir/transposonPSI.sh $BestAss
```  


** % bases maked by repeatmasker: 4.75%**

** % bases masked by transposon psi: 4.19% **


# Gene Prediction
Gene prediction followed two steps:
Pre-gene prediction - Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
Gene models were used to predict genes in the Fusarium genome. This used results from CEGMA as hints for gene models.

## Pre-gene prediction
Quality of genome assemblies was assessed by looking for the gene space in the assemblies.

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
  Assembly=repeat_masked/F.venenatum/strain1/filtered_contigs_repmask/strain1_contigs_unmasked.fa
  qsub $ProgDir/sub_cegma.sh $Assembly dna
```
The cegma completeness report gave an indication of the number of genes core
eukaryotic genes were present:
** Number of cegma genes present and complete: 237 (95.56%) **
** Number of cegma genes present and partial: 241 (97.18%) **

##Gene prediction

Gene prediction was performed for the neonectria genome.
CEGMA genes could be used as hints for the location of CDS.

For the moment we shall just use the gene model trained to Fusarium.
This model is from a closely related organism that is also plant pathogen.

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
  Assembly=repeat_masked/F.venenatum/strain1/filtered_contigs_repmask/strain1_contigs_unmasked.fa
  GeneModel=fusarium
  qsub $ProgDir/submit_augustus.sh $GeneModel $Assembly false
```

** Number of genes predicted: 11837 **


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

## Identification of potential CrisprCas sites.

The script optimus was used to identify potential CrisprCas sites Within
predicted proteins.

The commands to do this were:

```bash
  # GeneSeq=gene_pred/augustus/neonectria_galligena/NG-R0905_EMR/NG-R0905_EMR_aug_out.codingseq
  Organism=F.venenatum
  Strain=strain1
  ProgDir=~/git_repos/emr_repos/scripts/fusarium_venenatum/OPTIMus
  OutDir=analysis/protospacers/$Organism/$Strain
  GeneSeq=$(ls gene_pred/augustus/$Organism/$Strain/*_aug_out.codingseq)
  mkdir -p $OutDir
  $ProgDir/journal.pone.0133085.s004.pl $GeneSeq "threshold" 1 > $OutDir/"$Strain"_protospacer_sites.txt
  $ProgDir/Optimus2csv.py --inp $OutDir/"$Strain"_protospacer_sites.txt  --out $OutDir/"$Strain"_protospacer_by_gene.csv
```

It was realised that intron-exon boundaries would interfere with the prediction
of protospacers from cds. For this reason a new protospacer prediction program
was written.

```bash
  Organism=F.venenatum
  Strain=strain1
  ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
  Aug_Gff=gene_pred/augustus/$Organism/$Strain/"$Strain"_augustus_preds.gtf
  OutDir=analysis/protospacers/$Organism/$Strain
  AugDB=$OutDir/"$Strain"_Aug.db
  $ProgDir/make_gff_database.py --inp $Aug_Gff --db $AugDB
  ProgDir=~/git_repos/emr_repos/scripts/fusarium_venenatum/OPTIMus
  ProtospacerCSV=OutDir/"$Strain"_protospacer_by_gene.csv
  $ProgDir/protospacer_finder.py --inp $AugDB --out $OutDir/"$Strain"_protospacer_by_gene.csv
```
