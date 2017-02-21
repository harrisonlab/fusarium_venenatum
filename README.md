# fusarium_venenatum
Bioinformatic analysis of fusarium venenatum genomes

All work was carried out in the directory:

```bash
cd /home/groups/harrisonlab/project_files
mkdir -p fusarium_venenatum
ls /home/groups/harrisonlab/project_files/fusarium_venenatum
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


# Data organisation


## Building of directory structure
```bash
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium_venenatum
	mkdir -p $ProjectDir/raw_dna/paired/F.venenatum/strain1/F
	mkdir -p $ProjectDir/raw_dna/paired/F.venenatum/strain1/R
  RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160401_M004465_0007-AGKF2
  ProjectDir=/home/groups/harrisonlab/project_files/fusarium_venenatum
  mkdir -p $ProjectDir/raw_dna/paired/F.venenatum/C1/F
  mkdir -p $ProjectDir/raw_dna/paired/F.venenatum/C1/R
  mkdir -p $ProjectDir/raw_dna/paired/F.venenatum/C2/F
  mkdir -p $ProjectDir/raw_dna/paired/F.venenatum/C2/R
  mkdir -p $ProjectDir/raw_dna/paired/F.venenatum/C3/F
  mkdir -p $ProjectDir/raw_dna/paired/F.venenatum/C3/R
  mkdir -p $ProjectDir/raw_dna/paired/F.venenatum/C4/F
  mkdir -p $ProjectDir/raw_dna/paired/F.venenatum/C4/R
  mkdir -p $ProjectDir/raw_dna/paired/F.venenatum/C5/F
  mkdir -p $ProjectDir/raw_dna/paired/F.venenatum/C5/R
  mkdir -p $ProjectDir/raw_dna/paired/F.venenatum/C6/F
  mkdir -p $ProjectDir/raw_dna/paired/F.venenatum/C6/R
  mkdir -p $ProjectDir/raw_dna/paired/F.venenatum/WT/F
  mkdir -p $ProjectDir/raw_dna/paired/F.venenatum/WT/R
```
Sequence data was moved into the appropriate directories

```bash
RawDatDir=/home/miseq_data/.tmp_nas_data/miseq_data/miseq_data/RAW/2016/160304_M04465_0005_000000000-AKTC6/Data/Intensities/BaseCalls
ProjectDir=/home/groups/harrisonlab/project_files/fusarium_venenatum
cp $RawDatDir/C1_S2_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.venenatum/C1/F/.
cp $RawDatDir/C1_S2_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.venenatum/C1/R/.
cp $RawDatDir/C2_S3_L001_R1_001.fastq.gz  $ProjectDir/raw_dna/paired/F.venenatum/C2/F/.
cp $RawDatDir/C2_S3_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.venenatum/C2/R/.
cp $RawDatDir/C3_S4_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.venenatum/C3/F/.
cp $RawDatDir/C3_S4_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.venenatum/C3/R/.
cp $RawDatDir/C5_S5_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.venenatum/C5/F/.
cp $RawDatDir/C5_S5_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.venenatum/C5/R/.
cp $RawDatDir/FvenWT_S1_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.venenatum/WT/F/.
cp $RawDatDir/FvenWT_S1_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.venenatum/WT/R/.
  RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160401_M004465_0007-AGKF2
  ProjectDir=/home/groups/harrisonlab/project_files/fusarium_venenatum
	cp $RawDatDir/FvenC4_S3_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.venenatum/C4/F/.
	cp $RawDatDir/FvenC4_S3_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.venenatum/C4/R/.
  cp $RawDatDir/FvenC6_S4_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.venenatum/C6/F/.
  cp $RawDatDir/FvenC6_S4_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.venenatum/C6/R/.
  cp $RawDatDir/FvenWT_S2_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.venenatum/WT/F/.
  cp $RawDatDir/FvenWT_S2_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.venenatum/WT/R/.
  RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160415_M004465_00011-AMLCL
  ProjectDir=/home/groups/harrisonlab/project_files/fusarium_venenatum
  cp $RawDatDir/FvenWT_S3_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.venenatum/WT/F/.
  cp $RawDatDir/FvenWT_S3_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.venenatum/WT/R/.
```


#Data qc

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```bash
  for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz | grep -v 'strain1'); do
  	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
  	echo $RawData;
  	qsub $ProgDir/run_fastqc.sh $RawData
  done
```

Trimming was performed on data to trim adapters from
sequences and remove poor quality data. This was done with fastq-mcf

Firstly, those strains with more than one run were identified:

```bash
for Strain in $(ls -d raw_dna/paired/*/*); do
NumReads=$(ls $Strain/F/*.gz | wc -l);
if [ $NumReads -gt 1 ]; then
echo "$Strain";
echo "$NumReads";
fi;
done
```

```
raw_dna/paired/F.venenatum/strain1
2
raw_dna/paired/F.venenatum/WT
3
```

Trimming was first performed on all strains that had a single run of data:

```bash
for StrainPath in $(ls -d raw_dna/paired/*/* | grep -v -e 'strain1' -e 'WT'); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
ReadsF=$(ls $StrainPath/F/*.fastq*)
ReadsR=$(ls $StrainPath/R/*.fastq*)
echo $ReadsF
echo $ReadsR
qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
done
```


Trimming was then performed for strains with multiple runs of data

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
	IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
  echo "strain1"
  StrainPath=raw_dna/paired/F.venenatum/strain1
  ReadsF=$(ls $StrainPath/F/fungus1_S1_L001_R1_001.fastq.gz)
  ReadsR=$(ls $StrainPath/R/fungus1_S1_L001_R2_001.fastq.gz)
  qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	StrainPath=raw_dna/paired/F.venenatum/strain1
	ReadsF=$(ls $StrainPath/F/fungus2_S1_L001_R1_001.fastq.gz)
	ReadsR=$(ls $StrainPath/R/fungus2_S1_L001_R2_001.fastq.gz)
	qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	echo "WT"
  StrainPath=raw_dna/paired/F.venenatum/WT
  ReadsF=$(ls $StrainPath/F/FvenWT_S1_L001_R1_001.fastq.gz)
  ReadsR=$(ls $StrainPath/R/FvenWT_S1_L001_R2_001.fastq.gz)
  qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	StrainPath=raw_dna/paired/F.venenatum/WT
	ReadsF=$(ls $StrainPath/F/FvenWT_S2_L001_R1_001.fastq.gz)
	ReadsR=$(ls $StrainPath/R/FvenWT_S2_L001_R2_001.fastq.gz)
	qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	StrainPath=raw_dna/paired/F.venenatum/WT
	ReadsF=$(ls $StrainPath/F/FvenWT_S3_L001_R1_001.fastq.gz)
	ReadsR=$(ls $StrainPath/R/FvenWT_S3_L001_R2_001.fastq.gz)
	qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
```


Data quality was visualised once again following trimming:
```bash
  for RawData in $(ls qc_dna/paired/*/*/*/*.fq.gz); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    echo $RawData;
    qsub $ProgDir/run_fastqc.sh $RawData
  done
```

kmer counting was performed using kmc
This allowed estimation of sequencing depth and total genome size

This was performed for strains with single runs of data

```bash
	for TrimPath in $(ls -d qc_dna/paired/*/* | grep -v -e 'WT' -e 'strain1'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		TrimF=$(ls $TrimPath/F/*.fq.gz)
		TrimR=$(ls $TrimPath/R/*.fq.gz)
		echo $TrimF
		echo $TrimR
		qsub $ProgDir/kmc_kmer_counting.sh $TrimF $TrimR
	done
```

and for strains with multiple runs of data:

```bash
	for TrimPath in $(ls -d qc_dna/paired/*/strain1); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		TrimF1=$(ls $TrimPath/F/fungus1_S1_L001_R1_001.fq.gz)
		TrimR1=$(ls $TrimPath/R/fungus1_S1_L001_R2_001.fq.gz)
		echo $TrimF1
		echo $TrimR1
		TrimF2=$(ls $TrimPath/F/fungus2_S1_L001_R1_001.fq.gz)
		TrimR2=$(ls $TrimPath/R/fungus2_S1_L001_R2_001.fq.gz)
		echo $TrimF2
		echo $TrimR2
		qsub $ProgDir/kmc_kmer_counting.sh $TrimF1 $TrimR1 $TrimF2 $TrimR2
	done
	for TrimPath in $(ls -d qc_dna/paired/*/WT); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		TrimF1=$(ls $TrimPath/F/FvenWT_S2_L001_R1_001.fq.gz)
		TrimR1=$(ls $TrimPath/R/FvenWT_S2_L001_R2_001.fq.gz)
		echo $TrimF1
		echo $TrimR1
		TrimF2=$(ls $TrimPath/F/FvenWT_S3_L001_R1_001.fq.gz)
		TrimR2=$(ls $TrimPath/R/FvenWT_S3_L001_R2_001.fq.gz)
		echo $TrimF2
		echo $TrimR2
		qsub $ProgDir/kmc_kmer_counting.sh $TrimF1 $TrimR1 $TrimF2 $TrimR2
	done
```

mode kmer abundance prior to error correction was reported using the following
commands:

```bash
  for File in $(ls qc_dna/kmc/*/*/*_true_kmer_summary.txt); do
    basename $File;
    tail -n3 $File | head -n1 ;
  done
```

<!--
kmer counting was performed using kmc.
This allowed estimation of sequencing depth and total genome size:

```bash
  qsub /home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc/kmer_counting.sh qc_dna/paired/F.venenatum/strain1/F/strain1_qc_F.fastq qc_dna/paired/F.venenatum/strain1/R/strain1_qc_R.fastq qc_dna/paired/F.venenatum/strain1/kmer_count
```

** Estimated Genome Size is: **

** Esimated Coverage is: **


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
 -->


 #Assembly

 Assembly was performed with:
 * Spades

 ## Spades Assembly



 ```bash
 	for StrainPath in $(ls -d qc_dna/paired/*/* | grep -v -e 'strain1' -e 'WT'); do
 		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
 		Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
 		Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
 		F_Read=$(ls $StrainPath/F/*.fq.gz)
 		R_Read=$(ls $StrainPath/R/*.fq.gz)
 		OutDir=assembly/spades/$Organism/$Strain
 		Jobs=$(qstat | grep 'submit_SPA' | grep 'qw' | wc -l)
 		while [ $Jobs -gt 1 ]; do
 			sleep 5m
 			printf "."
 			Jobs=$(qstat | grep 'submit_SPA' | grep 'qw' | wc -l)
 		done		
 		printf "\n"
 		echo $F_Read
 		echo $R_Read
 		qsub $ProgDir/submit_SPAdes.sh $F_Read $R_Read $OutDir correct 20
 	done
 ```

 Assembly for strains failed due to a lack of memory, as such the assembly was
 resubmitted with more RAM.

 ```bash
 	for StrainPath in $(ls -d qc_dna/paired/*/* | grep -v -e 'strain1' -e 'WT'); do
 		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
 		Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
 		Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
 		F_Read=$(ls $StrainPath/F/*.fq.gz)
 		R_Read=$(ls $StrainPath/R/*.fq.gz)
 		OutDir=assembly/spades/$Organism/$Strain
 		echo $F_Read
 		echo $R_Read
 		qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $OutDir correct 20
 	done
 ```

 Assemblies were submitted for genomes with data from multiple sequencing runs:

 ```bash
  for StrainPath in $(ls -d qc_dna/paired/F.*/WT); do
    echo $StrainPath
    ProgDir=/home/ransoe/git_repos/tools/seq_tools/assemblers/spades/multiple_libraries
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    echo $Strain
    echo $Organism
    TrimF1=$(ls $TrimPath/F/FvenWT_S2_L001_R1_001.fastq.gz)
    TrimR1=$(ls $TrimPath/R/FvenWT_S2_L001_R2_001.fastq.gz)
    echo $TrimF1
    echo $TrimR1
    TrimF2=$(ls $TrimPath/F/FvenWT_S3_L001_R1_001.fastq.gz)
    TrimR2=$(ls $TrimPath/R/FvenWT_S3_L001_R2_001.fastq.gz)
    echo $TrimF2
    echo $TrimR2
    OutDir=assembly/spades/$Organism/$Strain
    qsub $ProgDir/subSpades_2lib.sh $TrimF1 $TrimR1 $TrimF2 $TrimR2 $OutDir correct 30
  done
 ```
<!--
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
``` -->

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

The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and hardmasked files.

```bash

for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```

# Gene Prediction

Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
	Gene model training
		- Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline
	Gene prediction
		- Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.


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


## Gene prediction 1 - Braker1 gene model training and prediction

Gene prediction was performed using Braker1.

First, RNAseq data was aligned to Fusarium genomes.
* Greg had aligned RNAseq data to the genome using his Tophat pipeline. The
Acceptedhits.bam files were used as evidence for gene models training using
Braker and CodingQuary.

#### Braker prediction

```bash
  for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    mkdir -p alignment/$Organism/$Strain/concatenated
    # samtools merge -f alignment/$Organism/$Strain/concatenated/concatenated.bam \
    # ../quorn/tophat/WTCHG_25*/accepted_hits.bam
    OutDir=gene_pred/braker/$Organism/"$Strain"_braker
    AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
    GeneModelName="$Organism"_"$Strain"_braker
    rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
    qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  done
```

** Number of genes predicted:  **


## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Firstly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
mkdir -p $OutDir
AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
done
```

Secondly, genes were predicted using CodingQuary:

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep -e 'strain1'); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		OutDir=gene_pred/codingquary/$Organism/$Strain
		CufflinksGTF=gene_pred/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
		qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
	done
```

Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

```bash
for BrakerGff in $(ls gene_pred/braker/F.*/*_braker/*/augustus.gff3); do
Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker//g')
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
# BrakerGff=gene_pred/braker/$Organism/$Strain/F.oxysporum_fsp_cepae_Fus2_braker/augustus_extracted.gff
Assembly=$(ls repeat_masked/$Organism/$Strain/*/"$Strain"_contigs_softmasked.fa)
CodingQuaryGff=gene_pred/codingquary/$Organism/$Strain/out/PredictedPass.gff3
PGNGff=gene_pred/codingquary/$Organism/$Strain/out/PGN_predictedPass.gff3
AddDir=gene_pred/codingquary/$Organism/$Strain/additional
FinalDir=gene_pred/final/$Organism/$Strain/final
AddGenesList=$AddDir/additional_genes.txt
AddGenesGff=$AddDir/additional_genes.gff
FinalGff=$AddDir/combined_genes.gff
mkdir -p $AddDir
mkdir -p $FinalDir

bedtools intersect -v -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
# GffFile=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/additional/additional_genes.gff
# GffFile=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/out/PredictedPass.gff3

$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary
cp $BrakerGff $FinalDir/final_genes_Braker.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta


GffBraker=$FinalDir/final_genes_CodingQuary.gff3
GffQuary=$FinalDir/final_genes_Braker.gff3
GffAppended=$FinalDir/final_genes_appended.gff3
cat $GffBraker $GffQuary > $GffAppended

# cat $BrakerGff $AddDir/additional_gene_parsed.gff3 | bedtools sort > $FinalGff
done
```

The final number of genes per isolate was observed using:
```bash
for DirPath in $(ls -d gene_pred/final/F.*/*/final); do
echo $DirPath;
cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
echo "";
done
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
screen -a
cd /home/groups/harrisonlab/project_files/fusarium_venenatum
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
for Genes in $(ls gene_pred/final/F.*/*/*/final_genes_combined.pep.fasta); do
echo $Genes
$ProgDir/sub_interproscan.sh $Genes
done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
  for Proteins in $(ls gene_pred/final/F.*/*/*/final_genes_combined.pep.fasta); do
    Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    echo $Strain
    InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
    $ProgDir/append_interpro.sh $Proteins $InterProRaw
  done
```

## B) SwissProt

```bash
  for Proteome in $(ls gene_pred/final/F.*/*/*/final_genes_combined.pep.fasta); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=gene_pred/swissprot/$Organism/$Strain
    SwissDbDir=../../uniprot/swissprot
    SwissDbName=uniprot_sprot
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
    qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
  done
```

```bash
  for SwissTable in $(ls gene_pred/swissprot/*/*/swissprot_v2015_10_hits.tbl); do
    Strain=$(echo $SwissTable | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $SwissTable | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutTable=gene_pred/swissprot/$Organism/$Strain/swissprot_vJul2016_tophit_parsed.tbl
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
    $ProgDir/swissprot_parser.py --blast_tbl $SwissTable --blast_db_fasta ../../uniprot/swissprot/uniprot_sprot.fasta > $OutTable
  done
```

#Genomic analysis
<!-- The first analysis was based upon BLAST searches for genes known to be involved in toxin production -->


## D) AntiSMASH

Antismash was run to identify clusters of secondary metabolite genes within
the genome. Antismash was run using the weserver at:
http://antismash.secondarymetabolites.org


Results of web-annotation of gene clusters within the assembly were downloaded to
the following directories:

```bash
  for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep 'strain1'); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=analysis/antismash/$Organism/$Strain
    mkdir -p $OutDir
  done
```

```bash
  for Zip in $(ls analysis/antismash/*/*/*.zip); do
    OutDir=$(dirname $Zip)
    unzip -d $OutDir $Zip
  done
```

```bash
  for AntiSmash in $(ls analysis/antismash/*/*/*/*.final.gbk); do
    Organism=$(echo $AntiSmash | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $AntiSmash | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/antismash/$Organism/$Strain
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
    $ProgDir/antismash2gff.py --inp_antismash $AntiSmash > $OutDir/"$Strain"_secondary_metabolite_regions.gff
    printf "Number of clusters detected:\t"
    cat $OutDir/"$Strain"_secondary_metabolite_regions.gff | grep 'antismash_cluster' | wc -l
    # GeneGff=gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3
    GeneGff=gene_pred/final/$Organism/$Strain/final/final_genes_appended.gff3
    bedtools intersect -u -a $GeneGff -b $OutDir/"$Strain"_secondary_metabolite_regions.gff > $OutDir/metabolite_cluster_genes.gff
    cat $OutDir/metabolite_cluster_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > $OutDir/metabolite_cluster_gene_headers.txt
    printf "Number of predicted proteins in clusters:\t"
    cat $OutDir/metabolite_cluster_gene_headers.txt | wc -l
    printf "Number of predicted genes in clusters:\t"
    cat $OutDir/metabolite_cluster_genes.gff | grep -w 'gene' | wc -l
  done
```

```
  F.venenatum - strain1
  Number of clusters detected:    38
  Number of predicted proteins in clusters:       585
  Number of predicted genes in clusters:  562
```

These clusters represented the following genes. Note that these numbers just
show the number of intersected genes with gff clusters and are not confirmed by
function



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
