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

### Identifing read depth

```bash
  for Reads in $(ls raw_dna/minion/*/*/*_pass.fastq.gz); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    qsub $ProgDir/sub_count_nuc.sh 38 $Reads
  done
  for Reads in $(ls qc_dna/paired/*/*/*/*_trim.fq.gz | grep 'WT'); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    qsub $ProgDir/sub_count_nuc.sh 38 $Reads
  done
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

```bash
  for Assembly in $(ls assembly/spades_minion/*/*/contigs.fasta); do
    echo "Filtering contigs smaller than 500bp"
    InDir=$(dirname $Assembly)
    OutDir=$InDir/filtered_contigs
    mkdir -p $OutDir
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/abyss
    $ProgDir/filter_abyss_contigs.py $Assembly 500 > $OutDir/contigs_min_500bp.fasta
  done
```



### Quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades_minion/*/*//filtered_contigs/*_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    echo "$Organism - $Strain"
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

<!--
## Nanopolish scaffolding

```bash
OutDir=assembly/nanopolish
mkdir -p $OutDir
  Reads=$(ls raw_dna/minion/*/*/*_pass.fastq.gz)
Spades=$(ls /home/groups/harrisonlab/project_files/fusarium_venenatum/assembly/spades/F.venenatum/WT/filtered_contigs/contigs_min_500bp.fasta)
bwa index $Spades
bwa mem -x ont2d -t 8 $Spades $Reads > $OutDir/aligned.bam
samtools sort $OutDir/aligned.bam -f $OutDir/reads.sorted.bam
samtools view -u $OutDir/aligned.bam | samtools sort - -f $OutDir/reads.sorted.bam
samtools index reads.sorted.bam

TMPDIR=/tmp/nanopolish
mkdir -p $TMPDIR
bwa index $TMPDIR/contigs_min_500bp.fasta
bwa mem -x ont2d -t 16 $TMPDIR/contigs_min_500bp.fasta $TMPDIR/WT_07-03-17_pass.fastq.gz | samtools sort - -o $TMPDIR/reads.sorted.bam


python nanopolish_makerange.py $TMPDIR/contigs_min_500bp.fasta | parallel --results nanopolish.results -P 8 \
    nanopolish variants --consensus polished.{1}.fa -w {1} -r $TMPDIR/WT_07-03-17_pass.fastq.gz -b $TMPDIR/reads.sorted.bam -g $TMPDIR/contigs_min_500bp.fasta -t 8 --min-candidate-frequency 0.1
```
 -->

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

checking using busco

```bash
#for Assembly in  $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa); do
for Assembly in $(ls assembly/merged_canu_spades/*/*_first/merged.fasta); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  # BuscoDB="Fungal"
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
  OutDir=gene_pred/busco/$Organism/$Strain/assembly
  qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```

```bash
  for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt); do  
    echo $File;
    cat $File | grep -e '(C)' -e 'Total';
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

Contigs were renamed in accordance with ncbi suggestions for exclusion of
contigs. The Cpontamination screen report was downloaded to location in NCBI_report
below and renamed to StrainName_ncbi_report.txt

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  for Assembly in $(ls assembly/merged_canu_spades/F.venenatum/WT_spades_first/polished/split_contigs/*_contigs_renamed.fasta); do
    Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev | cut -f1 -d '_')
    # NCBI_report=$(ls assembly/merged_canu_spades/$Organism/$Strain*/ncbi_report1/*report.txt)
    NCBI_report=report.txt
    # OutDir=$(dirname $NCBI_report)
    OutDir=assembly/merged_canu_spades/F.venenatum/WT_spades_first/polished/ncbi_report1
    mkdir -p $OutDir
    printf "Exclude:\nSequence name, length, apparent source\ncontig_8\t5513\tvector/etc" > $NCBI_report
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_ncbi_contigs_renamed.fasta --coord_file $NCBI_report
  done
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/merged_canu_spades/*/*/polished/ncbi_report1/*_contigs_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)  
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

# Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
for BestAss in $(ls assembly/merged_canu_spades/*/*/polished/ncbi_report1/*_contigs_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev | cut -f1 -d '_')
OutDir=repeat_masked/$Organism/"$Strain"_ncbi/ncbi_submission
qsub $ProgDir/rep_modeling.sh $BestAss $OutDir
qsub $ProgDir/transposonPSI.sh $BestAss $OutDir
done
```

```bash

for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep 'WT_ncbi'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa | grep 'WT_ncbi'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```

```bash
for RepDir in $(ls -d repeat_masked/F.*/*/*); do
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
<!--
```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
  Assembly=repeat_masked/F.venenatum/WT_ncbi/ncbi_submission/polished_contigs_softmasked_repeatmasker_TPSI_appended.fa
  qsub $ProgDir/sub_cegma.sh $Assembly dna
```
The cegma completeness report gave an indication of the number of genes core
eukaryotic genes were present:
** Number of cegma genes present and complete:  **
** Number of cegma genes present and partial:  **
-->

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep 'WT_ncbi'); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  # BuscoDB="Fungal"
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
* Greg had aligned RNAseq data to the genome using his STAR pipeline. The
Acceptedhits.bam files were used as evidence for gene models training using
Braker and CodingQuary.

Accepted hits .bam file were concatenated and indexed for use for gene model training:
<!--
```bash
ls /home/groups/harrisonlab/project_files/quorn/align/*.bam
BamFiles=$(ls /home/groups/harrisonlab/project_files/quorn/minion/*.Aligned.out.bam | grep -v 'SE' | tr -d '\n' | sed 's/.bam/.bam /g')
samtools merge -f alignment/$Organism/$Strain/concatenated/concatenated.bam $BamFiles
``` -->

```bash
OutDir=alignment/F.venenatum/WT/minion
mkdir -p $OutDir
for File in $(ls /home/groups/harrisonlab/project_files/quorn/minion/*.out.bam); do
Prefix=$(basename $File | sed 's/.bam//g')
echo $Prefix
samtools sort -o $File $Prefix > $OutDir/"$Prefix"_sorted.bam
done
BamFiles=$(ls $OutDir/*_sorted.bam | tr -d '\n' | sed 's/.bam/.bam /g')
samtools merge -f $OutDir/concatenated.bam $BamFiles
```

#### Braker prediction

```bash
  for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep 'WT_ncbi'); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    mkdir -p alignment/$Organism/$Strain/concatenated
    OutDir=gene_pred/braker/$Organism/"$Strain"_braker
    AcceptedHits=alignment/F.venenatum/WT/minion/concatenated.bam
    GeneModelName="$Organism"_"$Strain"_braker
    rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
    qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  done
```

** Number of genes predicted:  **


```bash
for BrakerGff in $(ls gene_pred/braker/F.*/*_braker/*/augustus.gff3  | grep 'WT_ncbi'); do
Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker//g')
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
FinalDir=gene_pred/final/$Organism/$Strain/final
Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_softmasked.fa)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
done

```

## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Firstly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep 'WT_ncbi'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
mkdir -p $OutDir
AcceptedHits=alignment/F.venenatum/WT/minion/concatenated.bam
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
done
```

Secondly, genes were predicted using CodingQuary:

```bash
  for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep 'WT_ncbi'); do
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
for BrakerGff in $(ls gene_pred/braker/F.*/*_braker/*/augustus.gff3  | grep 'WT_ncbi'); do
Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker//g')
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_softmasked.fa)
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
``` -->

## Assessing the Gene space in predicted transcriptomes:

```bash
	for Assembly in $(ls gene_pred/final/F.venenatum/WT_ncbi/final/final_genes_Braker.gene.fasta); do
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
  for File in $(ls gene_pred/busco/F*/*/genes/*/short_summary_*.txt); do  
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
  # for Genes in $(ls gene_pred/final/F.*/*/*/final_genes_combined.pep.fasta); do
  for Genes in $(ls gene_pred/braker/F.venenatum/WT_ncbi_braker/F.venenatum_WT_ncbi_braker/augustus.aa); do
    echo $Genes
    $ProgDir/sub_interproscan.sh $Genes
  done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
  # for Proteins in $(ls gene_pred/final/F.*/*/*/final_genes_combined.pep.fasta); do
  for Proteins in $(ls gene_pred/braker/F.venenatum/WT_ncbi_braker/F.venenatum_WT_ncbi_braker/augustus.aa); do
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
  # for Proteome in $(ls gene_pred/final/F.*/*/*/final_genes_combined.pep.fasta); do
  for Proteome in $(ls gene_pred/braker/F.venenatum/WT_ncbi_braker/F.venenatum_WT_ncbi_braker/augustus.aa); do
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
  for GeneGff in $(ls gene_pred/final/F.*/*/*/final_genes_appended.gff3); do
    Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
    Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_unmasked.fa)
    InterPro=$(ls gene_pred/interproscan/$Organism/$Strain/*_interproscan.tsv)
    SwissProt=$(ls gene_pred/swissprot/$Organism/$Strain/swissprot_vJul2016_tophit_parsed.tbl)
    OutDir=gene_pred/annotation/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/annotation_tables
    $ProgDir/build_annot_tab.py --genome $Assembly --genes_gff $GeneGff --InterPro $InterPro --Swissprot $SwissProt > $OutDir/"$Strain"_annotation.tsv
  done
``` -->


#Genomic analysis
<!-- The first analysis was based upon BLAST searches for genes known to be involved in toxin production -->


## D) Secondary metabolites (Antismash and SMURF)

Antismash was run to identify clusters of secondary metabolite genes within
the genome. Antismash was run using the weserver at:
http://antismash.secondarymetabolites.org

<!--
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

These clusters represented the following genes. Note that these numbers just
show the number of intersected genes with gff clusters and are not confirmed by
function

```
  F.venenatum - strain1
  Number of clusters detected:    38
  Number of predicted proteins in clusters:       585
  Number of predicted genes in clusters:  562
```
-->

SMURF was also run to identify secondary metabolite gene clusters.

Genes needed to be parsed into a specific tsv format prior to submission on the
SMURF webserver.

```bash
  # Gff=gene_pred/final/F.venenatum/strain1/final/final_genes_appended.gff3
  for Gff in $(ls gene_pred/final/F.venenatum/WT_ncbi/final/final_genes_Braker.gff3); do
    Strain=$(echo $Gff | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Gff | rev | cut -f4 -d '/' | rev)
    OutDir=analysis/secondary_metabolites/smurf/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
    $ProgDir/gff2smurf.py --gff $Gff > $OutDir/"$Strain"_genes_smurf.tsv
  done
```

SMURF output was received by email and downloaded to the cluster in the output
directory above.

Output files were parsed into gff format:

```bash
  for OutDir in $(ls -d analysis/secondary_metabolites/smurf/F.venenatum/WT_ncbi); do
    SmurfClusters=$(ls $OutDir/Secondary-Metabolite-Clusters.txt)
    SmurfBackbone=$(ls $OutDir/Backbone-genes.txt)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
    $ProgDir/smurf2gff.py --smurf_clusters $SmurfClusters --smurf_backbone $SmurfBackbone > $OutDir/Smurf_clusters.gff
  done
```
