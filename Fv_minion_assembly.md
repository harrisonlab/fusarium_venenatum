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

Checking assembly quality

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

## Identifying low coverage regions

Checking MiSeq coverage against WT contigs

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

Checking nanopore coverage against WT contigs

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
# Threshold=5
# FlaggedRegions=$AssemblyDir/aligned_MiSeq/pilon.fasta_flagged_regions_5x.txt
# ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
# $ProgDir/flag_low_coverage.py --genomecov $CoverageTxt --min $Threshold > $FlaggedRegions
```

Contigs with a low coverage from both nanopore and miseq reads:

```bash
  NanoporeFlagged=$(ls assembly/merged_canu_spades/F.venenatum/WT_spades_first/polished/aligned_minion/pilon.fasta_flagged_regions.txt)
  MiSeqFlagged=$(ls assembly/merged_canu_spades/F.venenatum/WT_spades_first/polished/aligned_MiSeq/pilon.fasta_flagged_regions.txt)
  ConsensusFlagged=
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
  $ProgDir/flag_low_cov_consensus.py --flagged1 $NanoporeFlagged --flagged2 $MiSeqFlagged > assembly/merged_canu_spades/F.venenatum/WT_spades_first/polished/consensus_flagged.txt
```

The output flagged regions were investigated by loading the assembly, Miseq
aligned reads and nanopore aligned reads into IGV. Each region was then visually
inspected to see the level of support from both alignments. Those regions without
coverage from either alignment were considered for splitting. The regions selected
for splitting(/removal) typically had long runs of a single or repetative nucleotides.
Results are saved in following excel spreadsheet and tab delimited file:

```bash
ls assembly/merged_canu_spades/F.venenatum/WT_spades_first/polished/consensus_flagged_region_edits.xlsx
ls assembly/merged_canu_spades/F.venenatum/WT_spades_first/polished/consensus_flagged_region_validity.txt
```

These files were parsed to make an instructions file for splitting contigs:

```bash
OutDir=assembly/merged_canu_spades/F.venenatum/WT_spades_first/polished/split_contigs
mkdir -p $OutDir
ConsensusValidity=assembly/merged_canu_spades/F.venenatum/WT_spades_first/polished/consensus_flagged_region_validity.txt
printf "Trim:\nSequence name, length, span(s), apparent source\n" > $OutDir/trim_instructions.txt
for Header in $(cat $ConsensusValidity | sed -e 's/\r/\n/g' | cut -f1,2 | grep -w 'y' | cut -f1 | cut -f1 -d ':' | sort | uniq); do
  ContigLgth=$(cat assembly/merged_canu_spades/F.venenatum/WT_spades_first/polished/consensus_flagged.txt | cut -f1,2 | grep -w "$Header" | cut -f2 | cut -f3 -d ' ')
  cat $ConsensusValidity | sed -e 's/\r/\n/g' | grep "$Header" | cut -f1,2 | grep -w 'y' | cut -f1 | cut -f2 -d ':' | grep '-' | sed 's/$/,/g' | tr -d '\n' | sed 's/-/../g' | sed -e "s/^/$Header\t$ContigLgth\t/g" | sed 's/,$/\treadmapping:no_coverage\n/g'
done >> $OutDir/trim_instructions.txt
```



## Editing contigs

Contigs were renamed in accordance with ncbi recomendations

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  for Assembly in $(ls assembly/merged_canu_spades/*/*/polished/pilon.fasta); do
    TrimInstructions=$(ls assembly/merged_canu_spades/F.venenatum/WT_spades_first/polished/split_contigs/trim_instructions.txt)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev | cut -f1 -d '_')
    OutDir=$(dirname $TrimInstructions)
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_contigs_renamed.fasta --coord_file $TrimInstructions
  done
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/merged_canu_spades/F.venenatum/WT_spades_first/polished/split_contigs/*_contigs_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

<!-- Contigs were renamed in accordance with ncbi suggestions for exclusion of
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
``` -->

# Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
for BestAss in $(ls assembly/merged_canu_spades/F.venenatum/WT_spades_first/polished/split_contigs/*_contigs_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev | cut -f1 -d '_')
OutDir=repeat_masked/$Organism/"$Strain"_ncbi/ncbi_submission
qsub $ProgDir/rep_modeling.sh $BestAss $OutDir
qsub $ProgDir/transposonPSI.sh $BestAss $OutDir
done
```
