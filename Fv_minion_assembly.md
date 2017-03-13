# Fusarium venenatum
==========

Scripts used for the analysis of Fusarium venenatum genomes
Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/fusarium_venenatum

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
  Organism=F.venenatum
  Strain=WT
  Date=07-03-17
  mkdir -p raw_dna/minion/$Organism/$Strain/$Date
  for Fast5Dir in $(ls -d $RawDatDir/*); do
    poretools fastq $Fast5Dir | gzip -cf
  done > raw_dna/minion/$Organism/$Strain/"$Strain"_"$Date"_pass.fastq.gz
  # poretools stats $RawDatDir/ > raw_dna/minion/$Organism/$Strain/"$Strain"_"$Date"_fail.stats.txt
  # poretools hist $RawDatDir/ > raw_dna/minion/$Organism/$Strain/"$Strain"_"$Date"_fail.hist
  # cat raw_dna/minion/$Organism/$Strain/"$Strain"_"$Date".fastq.gz raw_dna/minion/$Organism/$Strain/"$Strain"_"$Date"_fail.fastq.gz > raw_dna/minion/$Organism/$Strain/"$Strain"_"$Date"_pass-fail.fastq.gz
```


### Canu assembly

```bash
  Organism=F.venenatum
  Strain=WT
  Reads=$(ls raw_dna/minion/$Organism/$Strain/*_pass.fastq.gz)
  GenomeSz="38m"
  Prefix="$Strain"
  OutDir=assembly/canu-1.4/$Organism/"$Strain"
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  qsub $ProgDir/submit_canu.sh $Reads $GenomeSz $Prefix $OutDir
```


### Quast

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/canu-1.4/*/*/*.contigs.fasta); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
OutDir=assembly/canu-1.4/$Organism/$Strain/filtered_contigs
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```


<!-- Assemblies were polished using Pilon

```bash
for Assembly in $(ls assembly/canu-1.4/*/*/*.contigs.fasta); do
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
echo $Strain
echo $Organism
TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1 | tail -n1);
TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1 | tail -n1);
TrimF2_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n2 | tail -n1);
TrimR2_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n2 | tail -n1);
TrimF3_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n3 | tail -n1);
TrimR3_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n3 | tail -n1);
echo $TrimF1_Read
echo $TrimR1_Read
echo $TrimF2_Read
echo $TrimR2_Read
echo $TrimF3_Read
echo $TrimR3_Read
OutDir=assembly/canu-1.4/$Organism/$Strain/polished
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
qsub $ProgDir/sub_pilon_3_libs.sh $Assembly $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir
done
``` -->


This merged assembly was polished using Pilon

```bash
for Assembly in $(ls assembly/canu-1.4/*/*/*.contigs.fasta); do
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
echo $Strain
echo $Organism
TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n3 | tail -n1);
TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n3 | tail -n1);
  OutDir=assembly/canu-1.4/$Organism/$Strain/polished
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir
  done
```


<!--
Inspection of flagged regions didn't identify any contigs that needed to be broken.
-->


### Hybrid assembly:

#### Hybrid assembly: Spades Assembly

```bash
  for MinIonDat in $(ls raw_dna/minion/*/*/*_pass.fastq.gz); do
  Organism=$(echo $MinIonDat | rev | cut -f3 -d '/' | rev)
  Strain=$(echo $MinIonDat | rev | cut -f2 -d '/' | rev)
  IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
  echo $Strain
  echo $Organism
  TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1 | tail -n1);
  TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1 | tail -n1);
  TrimF2_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n2 | tail -n1);
  TrimR2_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n2 | tail -n1);
  TrimF3_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n3 | tail -n1);
  TrimR3_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n3 | tail -n1);
  echo $TrimF1_Read
  echo $TrimR1_Read
  echo $TrimF2_Read
  echo $TrimR2_Read
  echo $TrimF3_Read
  echo $TrimR3_Read
  OutDir=assembly/spades_minion/$Organism/$Strain
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades/multiple_libraries
  qsub $ProgDir/subSpades_3lib_minion.sh $MinIonDat $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir
  done
```


### Quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades_minion/*/*/contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    OutDir=assembly/spades_pacbio/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

# Merging Minion and Hybrid Assemblies

```bash
  for PacBioAssembly in $(ls assembly/canu-1.4/*/*/polished/*.fasta | grep -v -e '_nanopore' -e '_pass-fail'); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_*/$Organism/$Strain/contigs.fasta)
    OutDir=assembly/merged_canu_spades/$Organism/$Strain
    AnchorLength=500000
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
  done
```

```bash
  for PacBioAssembly in $(ls assembly/canu-1.4/*/*/polished/*.fasta | grep -v -e '_nanopore' -e '_pass-fail'); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_*/$Organism/$Strain/contigs.fasta)
    OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_100k_anchor
    AnchorLength=100000
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
  done
```

```bash
  for PacBioAssembly in $(ls assembly/canu-1.4/*/*/polished/*.fasta | grep -v -e '_nanopore' -e '_pass-fail'); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_*/$Organism/$Strain/contigs.fasta)
    OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_spades_first
    AnchorLength=500000
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $HybridAssembly $PacBioAssembly $OutDir $AnchorLength
  done
```

### Quast

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/merged_canu_spades/*/*/merged.fasta); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```


This merged assembly was polished using Pilon

```bash
for Assembly in $(ls assembly/merged_canu_spades/*/*/merged.fasta | grep 'spades_first'); do
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev | cut -f1 -d '_')
IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
echo $Strain
echo $Organism
TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n3 | tail -n1);
TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n3 | tail -n1);
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir/polished
done
```


Contigs were renamed in accordance with ncbi recomendations

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  # printf "contig_17\tsplit\t780978\t780971\tcanu:missassembly\n"
  for Assembly in $(ls assembly/merged_canu_spades/*/*/polished/pilon.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev | cut -f1 -d '_')
    OutDir=$(dirname $Assembly)
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_contigs_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

Contigs were renamed in accordance with ncbi suggestions for exclusion of
contigs. The Cpontamination screen report was downloaded to location in NCBI_report
below and renamed to StrainName_ncbi_report.txt

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  for Assembly in $(ls assembly/merged_canu_spades/*/*/polished/*_contigs_renamed.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev | cut -f1 -d '_')
    NCBI_report=$(ls assembly/merged_canu_spades/$Organism/$Strain*/ncbi_report1/*report.txt)
    OutDir=$(dirname $NCBI_report)
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_contigs_renamed.fasta --coord_file $NCBI_report
  done
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/merged_canu_spades/*/*/ncbi_report1/*_contigs_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

# Preliminary analysis

## Checking MiSeq coverage against WT contigs

```bash
for Assembly in $(ls assembly/merged_canu_spades/*/*/polished/pilon.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev | cut -f1 -d '_')
IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
echo $Strain
echo $Organism
TrimF1_Read=$(ls $IlluminaDir/F/*trim.fq.gz | head -n1 | tail -n1);
TrimR1_Read=$(ls $IlluminaDir/R/*trim.fq.gz | head -n1 | tail -n1);
TrimF2_Read=$(ls $IlluminaDir/F/*trim.fq.gz | head -n2 | tail -n1);
TrimR2_Read=$(ls $IlluminaDir/R/*trim.fq.gz | head -n2 | tail -n1);
TrimF3_Read=$(ls $IlluminaDir/F/*trim.fq.gz | head -n3 | tail -n1);
TrimR3_Read=$(ls $IlluminaDir/R/*trim.fq.gz | head -n3 | tail -n1);
echo $TrimF1_Read
echo $TrimR1_Read
echo $TrimF2_Read
echo $TrimR2_Read
echo $TrimF3_Read
echo $TrimR3_Read
# OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_414
InDir=$(dirname $Assembly)
OutDir=$InDir/aligned_MiSeq
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie_3lib.sh $Assembly $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir
done
```

## Identifying low coverage regions

 ```bash
  Assembly=$(ls assembly/merged_canu_spades/F.venenatum/WT_spades_first/polished/pilon.fasta)
  Reads=$(ls raw_dna/minion/*/*/*_pass.fastq.gz)
  AssemblyDir=$(dirname $Assembly)
  OutDir=$AssemblyDir/aligned_minion
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
  qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir
```

```bash
qlogin
cd /home/groups/harrisonlab/project_files/fusarium_venenatum
Assembly=$(ls assembly/merged_canu_spades/F.venenatum/WT_spades_first/polished/pilon.fasta)
AssemblyDir=$(dirname $Assembly)
AlignedBam=$(ls $AssemblyDir/aligned_MiSeq/pilon.fasta_aligned_sorted.bam)
CoverageTxt=$AssemblyDir/aligned_MiSeq/pilon.fasta_coverage.txt
bedtools genomecov -max 10 -d -ibam $AlignedBam -g $Assembly > $CoverageTxt
Threshold=10
FlaggedRegions=$AssemblyDir/aligned_MiSeq/pilon.fasta_flagged_regions.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
$ProgDir/flag_low_coverage.py --genomecov $CoverageTxt --min $Threshold > $FlaggedRegions
# Re-run with a lower Threshold
Threshold=5
FlaggedRegions=$AssemblyDir/aligned_MiSeq/pilon.fasta_flagged_regions_5x.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
$ProgDir/flag_low_coverage.py --genomecov $CoverageTxt --min $Threshold > $FlaggedRegions
```

# Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
for BestAss in $(ls assembly/merged_canu_spades/*/*/polished/*_contigs_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev | cut -f1 -d '_')
OutDir=repeat_masked/$Organism/"$Strain"_ncbi/ncbi_submission
qsub $ProgDir/rep_modeling.sh $BestAss $OutDir
qsub $ProgDir/transposonPSI.sh $BestAss $OutDir
done
```
