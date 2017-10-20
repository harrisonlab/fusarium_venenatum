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

MiSeq data was organised using:

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

Minion Seqeuncing data was extracted using commands documented in
Fv_minion_assembly.md, documented in this repository.



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


Find predicted coverage for these isolates:

```bash
for RawData in $(ls qc_dna/paired/*/*/*/*q.gz); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
GenomeSz=38
OutDir=$(dirname $RawData)
qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
done
```


```bash
  for StrainDir in $(ls -d qc_dna/paired/*/*); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls qc_dna/paired/*/"$Strain"/*/*.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```

```
C1	66.97
C2	45.69
C3	94.41
C4	125.6
C5	62
C6	79.56
WT	149.37
strain1	199.46
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

 Assemblies were submitted for genomes with data from multiple sequencing runs:

 ```bash
  for StrainPath in $(ls -d qc_dna/paired/F.*/strain1); do
    echo $StrainPath
    ProgDir=/home/ransoe/git_repos/tools/seq_tools/assemblers/spades/multiple_libraries
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    echo $Strain
    echo $Organism
    TrimF1=$(ls $StrainPath/F/*q.gz | head -n1 | tail -n1)
    TrimR1=$(ls $StrainPath/R/*q.gz | head -n1 | tail -n1)
    echo $TrimF1
    echo $TrimR1
    TrimF2=$(ls $StrainPath/F/*q.gz | head -n2 | tail -n1)
    TrimR2=$(ls $StrainPath/R/*q.gz | head -n2 | tail -n1)
    echo $TrimF2
    echo $TrimR2
    OutDir=assembly/spades/$Organism/"$Strain"_2
    qsub $ProgDir/subSpades_2lib.sh $TrimF1 $TrimR1 $TrimF2 $TrimR2 $OutDir correct
  done
 ```


 Quast

```bash
for Strain in $(ls -d assembly/spades/*/* | rev | cut -f1 -d'/' | rev); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
Assembly=$(ls -d assembly/spades/*/$Strain/filtered_contigs/contigs_min_500bp.fasta)
Species=$(echo $Assembly | rev | cut -f4 -d'/' | rev)
OutDir=$(ls -d assembly/spades/*/$Strain/filtered_contigs)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

The results of quast were shown using the following commands:

```bash
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/report.txt); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev);
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev);
    echo;
    echo $Organism;
    echo $Strain;
    cat $Assembly;
  done > assembly/quast_results.txt
```


A Bioproject and Biosample was made with NCBI genbank for submission of genomes.
Following the creation of these submissions, the .fasta assembly was uploaded
through the submission portal. A note was provided requesting that the assembly
be run through the contamination screen to aid a more detailed resubmission in
future. The returned FCSreport.txt was downloaded from the NCBI webportal and
used to correct the assembly to NCBI standards.

NCBI reports (FCSreport.txt) were manually downloaded to the following locations:

```bash
  for Assembly in $(ls assembly/spades/*/*/*/contigs_min_500bp.fasta | grep -v 'ncbi_edits' | grep -w 'WT'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    NCBI_report_dir=genome_submission/$Organism/$Strain/initial_submission
    mkdir -p $NCBI_report_dir
  done
```

These downloaded files were used to correct assemblies:

```bash
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep -w 'WT'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
NCBI_report=$(ls genome_submission/$Organism/$Strain/initial_submission/Contamination*.txt)
OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
mkdir -p $OutDir
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file $NCBI_report > $OutDir/log.txt
done
```

Quast

```bash
for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta | grep -w 'WT'); do
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta | grep -w 'WT'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```

```bash
  for File in $(ls gene_pred/assembly/F*/*/genes/*/short_summary_*.txt); do  
    echo $File;
    cat $File | grep -e '(C)' -e 'Total';
  done
```


# Repeatmasking

Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The best assembly was used to perform repeatmasking

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
  # BestAss=assembly/spades/F.venenatum/WT/filtered_contigs/contigs_min_500bp.fasta
  BestAss=assembly/spades/F.venenatum/WT/ncbi_edits/contigs_min_500bp_renamed.fasta
  OutDir=repeat_masked/F.venenatum/WT/illumina_assembly_ncbi
  qsub $ProgDir/rep_modeling.sh $BestAss $OutDir
  qsub $ProgDir/transposonPSI.sh $BestAss $OutDir
```  


** % bases maked by repeatmasker: 4.75%**

** % bases masked by transposon psi: 4.19% **

The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and hardmasked files.

```bash

for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep -w 'WT' | grep 'ncbi'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa | grep -w 'WT' | grep 'ncbi'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```

```bash
for RepDir in $(ls -d repeat_masked/F.*/*/* | grep -w 'WT' | grep 'ncbi'); do
Strain=$(echo $RepDir | rev | cut -f2 -d '/' | rev)
Organism=$(echo $RepDir | rev | cut -f3 -d '/' | rev)  
RepMaskGff=$(ls $RepDir/*_contigs_hardmasked.gff)
TransPSIGff=$(ls $RepDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
printf "$Organism\t$Strain\n"
# printf "The number of bases masked by RepeatMasker:\t"
sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
# printf "The number of bases masked by TransposonPSI:\t"
sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
# printf "The total number of masked bases are:\t"
cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
echo
done
```

F.venenatum	WT
302604
144657
438768

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

<!-- ```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
  Assembly=repeat_masked/F.venenatum/strain1/filtered_contigs_repmask/strain1_contigs_unmasked.fa
  qsub $ProgDir/sub_cegma.sh $Assembly dna
```
The cegma completeness report gave an indication of the number of genes core
eukaryotic genes were present:
** Number of cegma genes present and complete: 237 (95.56%) **
** Number of cegma genes present and partial: 241 (97.18%) ** -->

```bash
# for Assembly in $(ls assembly/spades/F.venenatum/WT/filtered_contigs/contigs_min_500bp.fasta); do
for Assembly in $(ls  repeat_masked/*/*/*/*_contigs_softmasked.fa | grep -w 'WT' | grep 'ncbi'); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
  OutDir=gene_pred/busco/$Organism/$Strain/assembly
  qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```

```bash
  for File in $(ls gene_pred/assembly/F*/*/genes/*/short_summary_*.txt); do  
    echo $File;
    cat $File | grep -e '(C)' -e 'Total';
  done
```



## Gene prediction 1 - Braker1 gene model training and prediction

Gene prediction was performed using Braker1.

First, RNAseq data was aligned to Fusarium genomes.
* Greg had aligned RNAseq data to the genome using his Tophat pipeline. The
Acceptedhits.bam files were used as evidence for gene models training using
Braker and CodingQuary.



#### Aligning

Insert sizes of the RNA seq library were unknown until a draft alignment could
be made. To do this tophat and cufflinks were run, aligning the reads against a
single genome. The fragment length and stdev were printed to stdout while
cufflinks was running.

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -w 'WT' | grep 'ncbi'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for FileF in $(ls ../quorn/filtered/*.1.fq | grep -v 'SE'); do
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
done
printf "\n"
FileR=$(echo $FileF | sed 's/.1.fq/.2.fq/g')
echo $FileF
echo $FileR
Prefix=$(echo $FileF | rev | cut -f1 -d '/' | rev | sed "s/.1.fq//g")
# Timepoint=$(echo $FileF | rev | cut -f2 -d '/' | rev)
Timepoint="treatment"
#echo "$Timepoint"
OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
done
done
```

Accepted hits .bam file were concatenated and indexed for use for gene model training:
<!--
```bash
ls /home/groups/harrisonlab/project_files/quorn/align/*.bam
BamFiles=$(ls /home/groups/harrisonlab/project_files/quorn/align/*.Aligned.sortedByCoord.out.bam | grep -v 'SE' | tr -d '\n' | sed 's/.bam/.bam /g')
samtools merge -f alignment/$Organism/$Strain/concatenated/concatenated.bam $BamFiles
``` -->

```bash
# For all alignments
BamFiles=$(ls alignment/star/F.venenatum/WT/treatment/*/*.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
OutDir=alignment/star/F.venenatum/WT/concatenated
mkdir -p $OutDir
samtools merge -f $OutDir/concatenated.bam $BamFiles
# one from each media type
BamFiles=$(ls alignment/star/F.venenatum/WT/treatment/WTCHG_25*_201/star_aligmentAligned.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
OutDir=alignment/star/F.venenatum/WT/concatenated
mkdir -p $OutDir
samtools merge -f $OutDir/one_per_media.bam $BamFiles

```

#### Braker prediction

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_softmasked_repeatmasker_TPSI_appended.fa | grep -w 'WT' | grep 'ncbi'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
mkdir -p alignment/$Organism/$Strain/concatenated
OutDir=gene_pred/braker/$Organism/"$Strain"_braker
AcceptedHits=$(ls alignment/star/F.venenatum/WT/concatenated/concatenated.bam)
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
  for Assembly in $(ls repeat_masked/*/*/*/*_softmasked_repeatmasker_TPSI_appended.fa | grep -w 'WT' | grep 'ncbi'); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
    # OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated2
    mkdir -p $OutDir
    AcceptedHits=$(ls alignment/star/F.venenatum/WT/concatenated/one_per_media.bam)
    # AcceptedHits=$(ls alignment/star/F.venenatum/WT/concatenated/concatenated.bam)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
  done
```

Secondly, genes were predicted using CodingQuary:

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_softmasked_repeatmasker_TPSI_appended.fa | grep -w 'WT' | grep 'ncbi'); do
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
for BrakerGff in $(ls gene_pred/braker/F.*/*_braker/*/augustus.gff3 | grep 'WT_braker'); do
Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker//g')
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
# BrakerGff=gene_pred/braker/$Organism/$Strain/F.oxysporum_fsp_cepae_Fus2_braker/augustus_extracted.gff
Assembly=$(ls repeat_masked/$Organism/*/*/*_softmasked_repeatmasker_TPSI_appended.fa | grep -w 'WT' | grep 'ncbi')
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


GffBraker=$FinalDir/final_genes_Braker.gff3
GffQuary=$FinalDir/final_genes_CodingQuary.gff3
GffAppended=$FinalDir/final_genes_appended.gff3
cat $GffBraker $GffQuary > $GffAppended

# cat $BrakerGff $AddDir/additional_gene_parsed.gff3 | bedtools sort > $FinalGff
done
```

In preperation for submission to ncbi, gene models were renamed and duplicate gene features were identified and removed.
 * no duplicate genes were identified


```bash
GffAppended=$(ls gene_pred/final/F.venenatum/WT/final/final_genes_appended.gff3)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/remove_dup_features.py --inp_gff $GffAppended

GffRenamed=$FinalDir/final_genes_appended_renamed.gff3
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/gff_rename_genes.py --inp_gff $GffAppended > $GffRenamed

Assembly=$(ls repeat_masked/$Organism/*/*/*_softmasked_repeatmasker_TPSI_appended.fa | grep -w 'WT' | grep 'ncbi')
$ProgDir/gff2fasta.pl $Assembly $GffRenamed gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed

# The proteins fasta file contains * instead of Xs for stop codons, these should
# be changed
sed -i 's/\*/X/g' gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.pep.fasta
```




The final number of genes per isolate was observed using:
```bash
for DirPath in $(ls -d gene_pred/final/F.*/*/final | grep -w 'WT'); do
echo $DirPath;
cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_appended_renamed.pep.fasta | grep '>' | wc -l;
echo "";
done
```

# Assessing gene space in predicted transcriptomes

```bash
for Assembly in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.gene.fasta | grep -w 'WT'); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  # BuscoDB="Fungal"
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
  OutDir=gene_pred/busco/$Organism/$Strain/genes
  qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```

```bash
  for File in $(ls gene_pred/busco/F*/*/genes/*/short_summary_*.txt | grep -w 'WT'); do  
    echo $File;
    cat $File | grep -e '(C)' -e 'Total';
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
for Genes in $(ls gene_pred/final/F.*/*/*/final_genes_appended_renamed.pep.fasta  | grep -w 'WT'); do
echo $Genes
$ProgDir/sub_interproscan.sh $Genes
done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
  for Proteins in $(ls gene_pred/final/F.*/*/*/final_genes_appended_renamed.pep.fasta | grep -w 'WT'); do
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
  for Proteome in $(ls gene_pred/final/F.*/*/*/final_genes_appended_renamed.pep.fasta | grep -w 'WT'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=gene_pred/swissprot/$Organism/$Strain
    SwissDbDir=../../uniprot/swissprot
    SwissDbName=uniprot_sprot
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
    qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
  done
```

<!--
## C) Summarising annotation in annotation table

```bash
  for GeneGff in $(ls gene_pred/final/F.*/*/*/final_genes_appended_renamed.gff3 | grep -w 'WT'); do
    Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
    Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_unmasked.fa | grep 'ncbi')
    InterPro=$(ls gene_pred/interproscan/$Organism/$Strain/*_interproscan.tsv)
    SwissProt=$(ls gene_pred/swissprot/$Organism/$Strain/swissprot_vJul2016_tophit_parsed.tbl)
    OutDir=gene_pred/annotation/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/annotation_tables
    $ProgDir/build_annot_tab.py --genome $Assembly --genes_gff $GeneGff --InterPro $InterPro --Swissprot $SwissProt > $OutDir/"$Strain"_annotation_ncbi.tsv
  done
```
 -->

#Genomic analysis
<!-- The first analysis was based upon BLAST searches for genes known to be involved in toxin production -->


## D) Secondary metabolites (Antismash and SMURF)

Antismash was run to identify clusters of secondary metabolite genes within
the genome. Antismash was run using the weserver at:
http://antismash.secondarymetabolites.org


Results of web-annotation of gene clusters within the assembly were downloaded to
the following directories:

```bash
  for Assembly in $(ls repeat_masked/*/*/illumina_assembly_ncbi/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w 'WT'); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
    mkdir -p $OutDir
  done
```

```bash
  for Zip in $(ls analysis/secondary_metabolites/antismash/*/*/*.zip); do
    OutDir=$(dirname $Zip)
    unzip -d $OutDir $Zip
  done
```

```bash
  for AntiSmash in $(ls analysis/secondary_metabolites/antismash/*/*/*/*.final.gbk); do
    Organism=$(echo $AntiSmash | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $AntiSmash | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
    Prefix=$OutDir/WT_antismash
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
    $ProgDir/antismash2gff.py --inp_antismash $AntiSmash --out_prefix $Prefix

    # Identify secondary metabolites within predicted clusters
    printf "Number of secondary metabolite detected:\t"
    cat "$Prefix"_secmet_clusters.gff | wc -l
    GeneGff=gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3
    bedtools intersect -u -a $GeneGff -b "$Prefix"_secmet_clusters.gff > "$Prefix"_secmet_genes.gff
    cat "$Prefix"_secmet_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > "$Prefix"_antismash_secmet_genes.txt
    bedtools intersect -wo -a $GeneGff -b "$Prefix"_secmet_clusters.gff | grep 'mRNA' | cut -f9,10,12,18 | sed "s/ID=//g" | perl -p -i -e "s/;Parent=g\w+//g" | perl -p -i -e "s/;Notes=.*//g" > "$Prefix"_secmet_genes.tsv
    printf "Number of predicted proteins in secondary metabolite clusters:\t"
    cat "$Prefix"_secmet_genes.txt | wc -l
    printf "Number of predicted genes in secondary metabolite clusters:\t"
    cat "$Prefix"_secmet_genes.gff | grep -w 'gene' | wc -l

      # Identify cluster finder additional non-secondary metabolite clusters
      printf "Number of cluster finder non-SecMet clusters detected:\t"
      cat "$Prefix"_clusterfinder_clusters.gff | wc -l
      GeneGff=gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3
      bedtools intersect -u -a $GeneGff -b "$Prefix"_clusterfinder_clusters.gff > "$Prefix"_clusterfinder_genes.gff
      cat "$Prefix"_clusterfinder_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > "$Prefix"_clusterfinder_genes.txt

      printf "Number of predicted proteins in cluster finder non-SecMet clusters:\t"
      cat "$Prefix"_clusterfinder_genes.txt | wc -l
      printf "Number of predicted genes in cluster finder non-SecMet clusters:\t"
      cat "$Prefix"_clusterfinder_genes.gff | grep -w 'gene' | wc -l
  done
```

These clusters represented the following genes. Note that these numbers just
show the number of intersected genes with gff clusters and are not confirmed by
function

```
F.venenatum - WT
Number of secondary metabolite detected:	35
Number of predicted proteins in secondary metabolite clusters:	986
Number of predicted genes in secondary metabolite clusters:	977
Number of cluster finder non-SecMet clusters detected:	86
Number of predicted proteins in cluster finder non-SecMet clusters:	2829
Number of predicted genes in cluster finder non-SecMet clusters:	2813
```

SMURF was also run to identify secondary metabolite gene clusters.

Genes needed to be parsed into a specific tsv format prior to submission on the
SMURF webserver.

```bash
  Gff=$(ls gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.gff3)
  OutDir=analysis/secondary_metabolites/smurf/F.venenatum/WT
  mkdir -p $OutDir
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
  $ProgDir/gff2smurf.py --gff $Gff > $OutDir/WT_genes_smurf.tsv
```

SMURF output was received by email and downloaded to the cluster in the output
directory above.

Output files were parsed into gff format:

```bash
  OutDir=analysis/secondary_metabolites/smurf/F.venenatum/WT
  Prefix="WT"
  GeneGff=gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.gff3
  SmurfClusters=$OutDir/Secondary-Metabolite-Clusters.txt
  SmurfBackbone=$OutDir/Backbone-genes.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
  $ProgDir/smurf2gff.py --smurf_clusters $SmurfClusters --smurf_backbone $SmurfBackbone > $OutDir/Smurf_clusters.gff
  bedtools intersect -wo -a $GeneGff -b $OutDir/Smurf_clusters.gff | grep 'mRNA' | cut -f9,10,12,18 | sed "s/ID=//g" | perl -p -i -e "s/;Parent=g\w+//g" | perl -p -i -e "s/;Notes=.*//g" > $OutDir/"$Prefix"_smurf_secmet_genes.tsv
```

Total number of secondary metabolite clusters:

```bash
for Assembly in $(ls repeat_masked/*/*/illumina_assembly_ncbi/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w 'WT'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
mkdir -p $OutDir
GeneGff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
AntismashClusters=$(ls analysis/secondary_metabolites/antismash/$Organism/$Strain/*_secmet_clusters.gff)
SmurfClusters=$(ls analysis/secondary_metabolites/smurf/$Organism/$Strain/Smurf_clusters.gff)
echo "Total number of Antismash clusters"
cat $AntismashClusters | wc -l
echo "Total number of SMURF clusters"
cat $SmurfClusters | wc -l
echo "number of Antismash clusters intersecting Smurf clusters"
bedtools intersect -a $AntismashClusters -b $SmurfClusters | wc -l
echo "number of Antismash clusters not intersecting Smurf clusters"
bedtools intersect -v -a $AntismashClusters -b $SmurfClusters | wc -l
echo "number of smurf clusters intersecting antismash clusters"
bedtools intersect -a $SmurfClusters -b $AntismashClusters | wc -l
echo "number of smurf clusters not intersecting antismash clusters"
bedtools intersect -v -a $SmurfClusters -b $AntismashClusters | wc -l
done
```


## E) Vitamin pathways

Of particular interest within primary metabolism are vitamin pathways.

Greg BLAST searched known Fg vitamin pathway genes (from KEGG analysis) against predicted gene models.

Vitamin pathway gene homolgs were provided in a table which was saved to the following file:

```bash
  mkdir analysis/vitamins
  less analysis/vitamins/Fg_vs_Fv_vitamin_hits.txt
  cat analysis/vitamins/Fg_vs_Fv_vitamin_hits.txt | sed 's/^M/newline/g' | sed 's/newline/\n/g' > analysis/vitamins/Fg_vs_Fv_vitamin_hits_parsed.txt
  cat analysis/vitamins/Fg_vs_Fv_vitamin_hits_parsed.txt | cut -f1 | grep -e "^g" > analysis/vitamins/Fv_vitamin_gene_headers.txt
```

These genes were extracted from annotation tables:

```bash
AnnotTab=$(ls gene_pred/annotation/F.venenatum/WT/WT_annotation_ncbi.tsv)
Vitamin_list=$(ls analysis/vitamins/Fv_vitamin_gene_headers.txt)
cat $Vitamin_list | sort | uniq | wc -l
cat $AnnotTab | grep -w -f $Vitamin_list > analysis/vitamins/Fv_vitamin_annotation.tsv
cat analysis/vitamins/Fv_vitamin_annotation.tsv | wc -l
```

# F) Genes with transcription factor annotations:


A list of PFAM domains, superfamily annotations used as part of the DBD database
and a further set of interproscan annotations listed by Shelest et al 2017 were made
http://www.transcriptionfactor.org/index.cgi?Domain+domain:all
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5415576/


<!-- ```bash
AnnotTab=$(ls gene_pred/annotation/F.venenatum/WT/WT_annotation_ncbi.tsv)
OutDir=analysis/transcription_factors
TF_Domains=$(ls /home/groups/harrisonlab/project_files/fusarium_venenatum/analysis/transcription_factors/TF_domains_DBD_Shelest_2017.txt)
cat $AnnotTab | grep -w -f $TF_Domains > $OutDir/TF_annotations.tsv
cat $OutDir/TF_annotations.tsv | cut -f1 > $OutDir/TF_headers.txt
cat $OutDir/TF_annotations.tsv | grep -e 'AS_Cluster' -e 'SM_Cluster' > $OutDir/TF_SecMet.tsv
``` -->
<!--
```bash
cat analysis/transcription_factors/DBD_SSF.tab | sed "s/\t/\tSSF/g" | awk '{ print $2 ", " $1}' | sed -r "s/^/(\'/g" | sed -r "s/$/\'),/g" | sed "s/, /\', \'/g" | sort | uniq | sed -e "s/ / /g" > analysis/transcription_factors/DBD.tsv
cat analysis/transcription_factors/DBD_PFAM.tab | awk '{ print $2 ", " $1}' | sed -r "s/^/(\'/g" | sed -r "s/$/\'),/g" | sed "s/, /\', \'/g" | sort | uniq | sed -e "s/ / /g" >> analysis/transcription_factors/DBD.tsv
cat analysis/transcription_factors/TFs_interpro_shelest_2017.txt | grep -v -e "^$" | sed -r "s/^/(\'/g" | sed -r "s/$/\'),/g" | sed "s/\t/\', \'/g" | sort | uniq | sed -e "s/ / /g" >> analysis/transcription_factors/DBD.tsv
``` -->

```bash
  for Interpro in $(ls gene_pred/interproscan/*/*/*_interproscan.tsv | grep -w 'WT'); do
    Organism=$(echo $Interpro | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Interpro | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/transcription_factors/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transcription_factors
    $ProgDir/interpro2TFs.py --InterPro $Interpro > $OutDir/"$Strain"_TF_domains.tsv
    echo "total number of transcription factors"
    cat $OutDir/"$Strain"_TF_domains.tsv | cut -f1 | sort | uniq > $OutDir/"$Strain"_TF_gene_headers.txt
    cat $OutDir/"$Strain"_TF_gene_headers.txt | wc -l
  done
```


## G) Summarising annotation in annotation table

```bash
for GeneGff in $(ls gene_pred/final/F.*/*/*/final_genes_appended_renamed.gff3 | grep -w 'WT'); do
Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_unmasked.fa | grep 'ncbi')
Antismash=$(ls analysis/secondary_metabolites/antismash/F.venenatum/WT/WT_antismash_secmet_genes.tsv)
Smurf=$(ls analysis/secondary_metabolites/smurf/F.venenatum/WT/WT_smurf_secmet_genes.tsv)
Vitamins=$(ls analysis/vitamins/Fg_vs_Fv_vitamin_hits_parsed.txt)
TFs=$(ls analysis/transcription_factors/F.venenatum/WT/WT_TF_domains.tsv)
InterPro=$(ls gene_pred/interproscan/$Organism/$Strain/*_interproscan.tsv)
SwissProt=$(ls gene_pred/swissprot/$Organism/$Strain/swissprot_vJul2016_tophit_parsed.tbl)
PH1_orthology=$(ls analysis/orthology/orthomcl/Fv_vs_Fg/Fv_vs_Fg_orthogroups.txt)
GR1_orthology=$(ls analysis/orthology/orthomcl/Fv_vs_Fg_JGI/Fv_vs_Fg_JGI_orthogroups.txt)
OutDir=gene_pred/annotation/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_venenatum/analysis/annotation_tables
$ProgDir/build_annot_Fv.py --genome $Assembly --genes_gff $GeneGff --Antismash $Antismash --Smurf $Smurf --vitamins $Vitamins --TFs $TFs --InterPro $InterPro --Swissprot $SwissProt --orthogroups_PH1 $PH1_orthology --orthogroups_GR1 $GR1_orthology > $OutDir/"$Strain"_annotation_ncbi.tsv
done
```

```bash
for AnnotTab in $(ls gene_pred/annotation/F.venenatum/WT/WT_annotation_ncbi.tsv); do
  Strain=$(echo $AnnotTab | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $AnnotTab | rev | cut -f3 -d '/' | rev)
  OutDir=$(ls -d analysis/transcription_factors/$Organism/$Strain)
  echo "Total number of clusters:"
  cat $AnnotTab | grep 'SecMet_cluster_' | cut -f7 | sort | uniq | wc -l
  cat $AnnotTab | grep 'SecMet_cluster_' | cut -f1,7 | cut -f1 > $OutDir/WT_SecMet_headers.txt
  cat $OutDir/WT_SecMet_headers.txt | cut -f1 -d '.' | sort | uniq | wc -l
  cat $AnnotTab | grep -e 'AS_' -e 'SM_' | cut -f1,14 | grep -v -e "\s$" | cut -f1 > $OutDir/WT_TF_SecMet_headers.txt
  cat $OutDir/WT_TF_SecMet_headers.txt | wc -l
done
```


<!--
##Genes with homology to PHIbase
Predicted gene models were searched against the PHIbase database using tBLASTx.

```bash

```

Top BLAST hits were used to annotate gene models.

```bash

```

following blasting PHIbase to the genome, the hits were filtered by effect on
virulence. -->

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

# Primary metabolism

Genes involved in primary metabolism can be identified through GO annotations:
http://amigo.geneontology.org/goose?query=SELECT+DISTINCT+descendant.acc%2C+descendant.name%2C+descendant.term_type%0D%0AFROM%0D%0A+term%0D%0A+INNER+JOIN+graph_path+ON+%28term.id%3Dgraph_path.term1_id%29%0D%0A+INNER+JOIN+term+AS+descendant+ON+%28descendant.id%3Dgraph_path.term2_id%29%0D%0AWHERE+term.name%3D%27Primary+Metabolic+Process%27+AND+distance+%3C%3E+0%3B&mirror=ebi&limit=0

GO terms could be searched in annotation files to look at these genes.
