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

  # Oxford nanopore 18/07/17 and other runs during HortRes conference
  RawDatDir=/home/miseq_data/minion/2017/*_FvenenatumWT/fast5/pass
  Organism=F.venenatum
  Strain=WT
  Date=18-07-17
  mkdir -p raw_dna/minion/$Organism/$Strain/$Date
  for Fast5Dir in $(ls -d $RawDatDir/*); do
    poretools fastq $Fast5Dir | gzip -cf
  done > raw_dna/minion/$Organism/$Strain/"$Strain"_"$Date"_pass.fastq.gz
```

Data was basecalled again using Albacore 2.02 on the minion server:

```bash
Organism=F.venenatum
Strain=WT
OutDir=/home/groups/harrisonlab/project_files/fusarium_venenatum/raw_dna/minion/$Organism/$Strain
mkdir -p $OutDir

ssh nanopore@nanopore
mkdir Fven_26-10-17
cd Fven_26-10-17
screen -a
# Oxford nanopore 07/03/17
Organism=F.venenatum
Strain=WT
Date=07-03-17
FlowCell="FLO-MIN106"
Kit="SQK-LSK108"
RawDatDir=/data/seq_data/minion/2017/Fvenenatum/downloaded/pass
# OutDir=/data/seq_data/minion/2017/Fvenenatum/albacore_v2.02
OutDir=/home/groups/harrisonlab/project_files/fusarium_venenatum/raw_dna/minion/$Organism/$Strain

mkdir -p ~/Fven_26-10-17/$Date
cd ~/Fven_26-10-17/$Date
~/.local/bin/read_fast5_basecaller.py \
  --flowcell $FlowCell \
  --kit $Kit \
  --input $RawDatDir \
  --recursive \
  --worker_threads 24 \
  --save_path "$Organism"_"$Strain"_"$Date" \
  --output_format fastq,fast5 \
  --reads_per_fastq_batch 4000
  cat "$Organism"_"$Strain"_"$Date"/workspace/pass/*.fastq | gzip -cf > "$Organism"_"$Strain"_"$Date"_albacore_v2.02.fastq.gz
  scp "$Organism"_"$Strain"_"$Date"_albacore_v2.02.fastq.gz armita@192.168.1.200:$OutDir/.
  tar -cz -f "$Organism"_"$Strain"_"$Date".tar.gz "$Organism"_"$Strain"_"$Date"
  mkdir -p $OutDir
  scp "$Organism"_"$Strain"_"$Date".tar.gz armita@192.168.1.200:$OutDir/.

  # Oxford nanopore 18/07/17
  RawDatDir=/data/seq_data/minion/2017/20170718_1517_FvenenatumWT/fast5/pass
  Organism=F.venenatum
  Strain=WT
  Date=18-07-17
  FlowCell="FLO-MIN107"
  Kit="SQK-LSK108"
  RawDatDir=/data/seq_data/minion/2017/20170718_1517_FvenenatumWT/fast5/pass
  OutDir=/home/groups/harrisonlab/project_files/fusarium_venenatum/raw_dna/minion/$Organism/$Strain

  mkdir -p ~/Fven_26-10-17/$Date
  cd ~/Fven_26-10-17/$Date
  ~/.local/bin/read_fast5_basecaller.py \
    --flowcell $FlowCell \
    --kit $Kit \
    --input $RawDatDir \
    --recursive \
    --worker_threads 12 \
    --save_path "$Organism"_"$Strain"_"$Date" \
    --output_format fastq,fast5 \
    --reads_per_fastq_batch 4000
    cat "$Organism"_"$Strain"_"$Date"/workspace/pass/*.fastq | gzip -cf > "$Organism"_"$Strain"_"$Date"_albacore_v2.02.fastq.gz
    scp "$Organism"_"$Strain"_"$Date"_albacore_v2.02.fastq.gz armita@192.168.1.200:$OutDir/.
  tar -cz -f "$Organism"_"$Strain"_"$Date".tar.gz "$Organism"_"$Strain"_"$Date"
  mkdir -p $OutDir
  mv "$Organism"_"$Strain"_"$Date".tar.gz $OutDir/.
  scp "$Organism"_"$Strain"_"$Date".tar.gz armita@192.168.1.200:$OutDir/.
```


### Identifing read depth

```bash
  for Reads in $(ls raw_dna/minion/*/*/*.fastq.gz | grep 'albacore' | grep '18-07-17'); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    qsub $ProgDir/sub_count_nuc.sh 38 $Reads
  done
  for Reads in $(ls raw_dna/paired/*/*/*/*.fq.gz | grep 'WT'); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    qsub $ProgDir/sub_count_nuc.sh 38 $Reads
  done
```


Splitting reads and trimming adapters using porechop
```bash
for RawReads in $(ls raw_dna/minion/*/*/*.fastq.gz | grep 'albacore' | grep '18-07-17'); do
Strain=$(echo $RawReads | rev | cut -f2 -d '/' | rev)
Organism=$(echo $RawReads | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=qc_dna/minion/$Organism/$Strain
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
qsub $ProgDir/sub_porechop.sh $RawReads $OutDir
done
```

Read coverage was estimated from the trimmed datasets:

```bash
GenomeSz=38
for Reads in $(ls qc_dna/minion/*/*/*.fastq.gz | grep 'albacore_v2.02'); do
echo $Reads
OutDir=$(dirname $Reads)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
qsub $ProgDir/sub_count_nuc.sh $GenomeSz $Reads $OutDir
done
for Reads in $(ls qc_dna/paired/*/*/*/*_trim.fq.gz | grep 'WT'); do
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
  qsub $ProgDir/sub_count_nuc.sh 38 $Reads
done
```

```bash
  for StrainDir in $(ls -d qc_dna/minion/F.venenatum/WT/); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls $StrainDir/*cov.txt | grep -v 'appended'); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```

```
WT	45.04
```

### Canu assembly

```bash
  Organism=F.venenatum
  Strain=WT
  Jobs=$(qstat | grep 'sub_porech' | grep 'qw' | wc -l)
  while [ $Jobs -gt 0 ]; do
  sleep 1m
  printf "."
  Jobs=$(qstat | grep 'sub_porech' | grep 'qw' | wc -l)
  done		
  printf "\n"
  # Reads=$(ls raw_dna/minion/$Organism/$Strain/*_pass.fastq.gz)
  Reads1=$(ls qc_dna/minion/$Organism/$Strain/*_trim.fastq.gz | grep 'albacore_v2.02' | head -n1 | tail -n1)
  Reads2=$(ls qc_dna/minion/$Organism/$Strain/*_trim.fastq.gz | grep 'albacore_v2.02' | head -n2 | tail -n1)
  GenomeSz="38m"
  Prefix="$Strain"
  OutDir=assembly/canu-1.6/$Organism/"$Strain"_albacore_v2
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  qsub $ProgDir/submit_canu_minion_2lib.sh $Reads1 $Reads2 $GenomeSz $Prefix $OutDir
```



### Assembbly using SMARTdenovo

```bash
for CorrectedReads in $(ls assembly/canu-1.6/F.venenatum/WT_albacore_v2/WT.trimmedReads.fasta.gz); do
Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
Prefix="$Strain"
OutDir=assembly/SMARTdenovo/$Organism/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/SMARTdenovo
qsub $ProgDir/sub_SMARTdenovo.sh $CorrectedReads $Prefix $OutDir
done
```

Contigs shorter than 500bp were removed from the assembly

```bash
for Contigs in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg | grep 'albacore'); do
AssemblyDir=$(dirname $Contigs)
mkdir $AssemblyDir/filtered_contigs
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/abyss
$ProgDir/filter_abyss_contigs.py $Contigs 500 > $AssemblyDir/contigs_min_500bp.fasta
done
```

<!-- Just running canu read correction:

```bash
  Organism=F.venenatum
  Strain=WT
  Reads=$(ls raw_dna/minion/$Organism/$Strain/*_pass.fastq.gz)
  GenomeSz="38m"
  Prefix="$Strain"
  OutDir=assembly/canu-1.6/$Organism/"$Strain"_correction
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  qsub $ProgDir/sub_canu_correction.sh $Reads $GenomeSz $Prefix $OutDir
``` -->
<!--
Running PBcR read correction (basis of canu read correcttion) using illumina data

```bash
qlogin -pe smp 8
CurDir=/home/groups/harrisonlab/project_files/fusarium_venenatum
cd $CurDir

Organism=F.venenatum
Strain=WT
Reads=$(ls $CurDir/raw_dna/minion/$Organism/$Strain/*_pass.fastq.gz)
IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
#TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1 | tail -n1);
#TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1 | tail -n1);
TrimF2_Read=$(ls $CurDir/$IlluminaDir/F/*_trim.fq.gz | head -n2 | tail -n1);
TrimR2_Read=$(ls $CurDir/$IlluminaDir/R/*_trim.fq.gz | head -n2 | tail -n1);
TrimF3_Read=$(ls $CurDir/$IlluminaDir/F/*_trim.fq.gz | head -n3 | tail -n1);
TrimR3_Read=$(ls $CurDir/$IlluminaDir/R/*_trim.fq.gz | head -n3 | tail -n1);
GenomeSz="38000000"
Prefix="$Strain"
OutDir=assembly/canu-1.6/$Organism/"$Strain"_correction

WorkDir=/tmp/PBcR
mkdir $WorkDir
cd $WorkDir

cp $Reads readsPacBio.fq.gz
gunzip readsPacBio.fq.gz
cp readsPacBio.fq | cut -f1 -d ' ' > readsPacBio2.fq

cp $TrimF3_Read readsF.fq.gz
cp $TrimR3_Read readsR.fq.gz
gunzip readsF.fq.gz
gunzip readsR.fq.gz

printf "merSize=14\n" > $Prefix.spec

fastqToCA -libraryname illumina -technology illumina-long -innie -insertsize 400 150 -mates $WorkDir/readsF.fq,$WorkDir/readsR.fq > $WorkDir/illumina.frg
# fastqToCA -libraryname illumina -technology illumina-long -innie -reads $WorkDir/readsF.fq  > $WorkDir/illumina.frg
PBcR -t 8 -length 500 -partitions 200 -l $Prefix -s $Prefix.spec genomeSize=$GenomeSz -fastq readsPacBio2.fq illumina.frg 2>&1 | tee log.txt

```
-->


### Quast

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
# for Assembly in $(ls assembly/canu-1.6/*/*/*.contigs.fasta | grep 'albacore'); do
for Assembly in $(ls assembly/SMARTdenovo/*/*/contigs_min_500bp.fasta | grep 'albacore'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
# OutDir=assembly/canu-1.6/$Organism/$Strain/filtered_contigs
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

Error correction using racon:

```bash
# for Assembly in $(ls assembly/canu-1.6/*/*/*.contigs.fasta | grep 'albacore'); do
for Assembly in $(ls assembly/SMARTdenovo/*/*/contigs_min_500bp.fasta | grep 'albacore'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev | sed 's/_albacore_v2//g')
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ReadsFq1=$(ls qc_dna/minion/$Organism/$Strain/*_trim.fastq.gz | grep 'albacore_v2.02' | head -n1 | tail -n1)
ReadsFq2=$(ls qc_dna/minion/$Organism/$Strain/*_trim.fastq.gz | grep 'albacore_v2.02' | head -n2 | tail -n1)
ReadsAppended=qc_dna/minion/$Organism/$Strain/"$Strain"_reads_appended.fastq.gz
cat $ReadsFq1 $ReadsFq2 > $ReadsAppended
# OutDir=assembly/canu-1.6/$Organism/$Strain/racon
OutDir=$(dirname $Assembly)/racon
Iterations=10
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/racon
# qsub $ProgDir/sub_racon.sh $Assembly $ReadsAppended $Iterations $OutDir
qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq1 $Iterations $OutDir
done
# rm $ReadsAppended
```

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/SMARTdenovo/*/*/racon/contigs_min_500bp_racon_round_10.fasta | grep 'WT' | grep 'round_10' | grep 'albacore'); do
    OutDir=$(dirname $Assembly)
    echo "" > tmp.txt
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    $ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/racon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
  done
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon/contigs_min_500bp_racon_round_10.fasta | grep 'WT' | grep 'round_10' | grep 'albacore'); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```


```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon/contigs_min_500bp_racon_round_10.fasta | grep 'WT' | grep 'round_10' | grep 'albacore'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt | grep 'WT'); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

# Assembly correction using nanopolish

Fast5 files are very large and need to be stored as gzipped tarballs. These needed temporarily unpacking but must be deleted after nanpolish has finished running.
<!--
The minion device lost connection to metrichor during it's run before reconnecting,
this is though to have led to the duplication of four reads in the dataset,
which prevented nanopolish from running. As such, these were removed during the
nanopolish extract step.
reads were:
```
>1c8314d1-022d-4745-8c82-f3ba3d4deefa_Basecall_Alignment_template:1D_000:template
>ab4ed1d5-a7b5-4d8b-ab2b-c0573f48be7d_Basecall_Alignment_template:1D_000:template
>bdfcaf0b-fc25-4413-ab79-46e6e99c9e1c_Basecall_Alignment_template:1D_000:template
>d6d26487-5839-4f57-908c-f170cc713971_Basecall_Alignment_template:1D_000:template
``` -->

Raw reads were moved onto the cluster scratch space for this step and unpacked:

```bash
ScratchDir=/data/scratch/nanopore_tmp_data/Fven
mkdir -p $ScratchDir
cp raw_dna/minion/F.venenatum/WT/*.tar.gz $ScratchDir/.
for Tar in $(ls $ScratchDir/*.tar.gz); do
  tar -zxvf $Tar -C $ScratchDir
done
```


```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon/racon_min_500bp_renamed.fasta | grep 'WT' | grep 'albacore'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Note - the full path from home must be used
ReadDir=raw_dna/nanopolish/$Organism/$Strain
# if [ -d $ReadDir ]; then
# echo "reads already extracted"
# else
# echo "extracting reads"
mkdir -p $ReadDir
# CurDir=$PWD
# cd $ReadDir
# Event information would have been used from all of the runs, howver MinKnow doesnt
# produce event-level information and therefore just the albacore data was used.
# for Fast5Dir in $(ls -d /home/miseq_data/minion/2017/Fvenenatum/downloaded/pass); do
# nanopolish extract -r $Fast5Dir \
# | gzip -cf
# done > "$Strain"_reads.fa.gz
# cd $CurDir
ReadsFq1=$(ls raw_dna/minion/F.venenatum/WT/F.venenatum_WT_07-03-17_albacore_v2.02.fastq.gz)
ReadsFq2=$(ls raw_dna/minion/F.venenatum/WT/F.venenatum_WT_18-07-17_albacore_v2.02.fastq.gz)
cat $ReadsFq1 $ReadsFq2 | gunzip -cf > $ReadDir/"$Strain"_concatenated_reads.fastq
/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish/nanopolish_remove_dup_reads.py --fastq $ReadDir/"$Strain"_concatenated_reads.fastq --out $ReadDir/"$Strain"_concatenated_reads_filtered.fastq

# cat $ReadsFq1 | gunzip -cf > $ReadDir/"$Strain"_07-03-17_reads.fastq
# cat $ReadsFq2 | gunzip -cf > $ReadDir/"$Strain"_18-07-17_reads.fastq
# /home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish/nanopolish_remove_dup_reads.py --fastq $ReadDir/"$Strain"_07-03-17_reads.fastq --out $ReadDir/"$Strain"_07-03-17_reads_filtered.fastq
# /home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish/nanopolish_remove_dup_reads.py --fastq $ReadDir/"$Strain"_18-07-17_reads.fastq --out $ReadDir/"$Strain"_18-07-17_reads_filtered.fastq

ScratchDir=/data/scratch/nanopore_tmp_data/Fven
Fast5Dir1=$ScratchDir/F.venenatum_WT_07-03-17/workspace/pass
Fast5Dir2=$ScratchDir/F.venenatum_WT_18-07-17/workspace/pass
# nanopolish extract -r -q -o $ReadDir/F.venenatum_WT_07-03-17_albacore_v2.02.fastq $Fast5Dir1
nanopolish index -d $Fast5Dir1 -d $Fast5Dir2 $ReadDir/"$Strain"_concatenated_reads_filtered.fastq

# $ReadDir/"$Strain"_concatenated_reads.fastq.gz | gunzip -cf \
# | grep -A3 \
# -e '1c8314d1-022d-4745-8c82-f3ba3d4deefa_Basecall_Alignment_template' \
# -e 'ab4ed1d5-a7b5-4d8b-ab2b-c0573f48be7d_Basecall_Alignment_template' \
# -e 'bdfcaf0b-fc25-4413-ab79-46e6e99c9e1c_Basecall_Alignment_template' \
# -e 'd6d26487-5839-4f57-908c-f170cc713971_Basecall_Alignment_template' \
# > tmp.txt
# cat $ReadDir/"$Strain"_reads.fa.gz | gunzip -cf \
# | grep -v -f tmp.txt | gzip -cf \
# > $ReadDir/"$Strain"_nanopolish_index.fastq.gz
# fi


# RawReads=$(ls $ReadDir/"$Strain"_reads_no_duplicates.fa.gz)
OutDir=$(dirname $Assembly)
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
qsub $ProgDir/sub_bwa_nanopolish.sh $Assembly $ReadDir/"$Strain"_concatenated_reads_filtered.fastq.fa.gz $OutDir/nanopolish
done
```

 Split the assembly into 50Kb fragments an submit each to the cluster for
 nanopolish correction

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon/racon_min_500bp_renamed.fasta | grep 'WT' | grep 'albacore'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)
RawReads=$(ls raw_dna/nanopolish/$Organism/$Strain/"$Strain"_concatenated_reads_filtered.fastq)
AlignedReads=$(ls $OutDir/nanopolish/reads.sorted.bam)

NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_makerange.py $Assembly > $OutDir/nanopolish/nanopolish_range.txt

Ploidy=1
echo "nanopolish log:" > nanopolish_log.txt
for Region in $(cat $OutDir/nanopolish/nanopolish_range.txt | tail -n+2); do
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
done		
printf "\n"
echo $Region
echo $Region >> nanopolish_log.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
qsub $ProgDir/sub_nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $Region $OutDir/$Region
done
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon/racon_min_500bp_renamed.fasta | grep 'WT' | grep 'albacore'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
OutDir=assembly/SMARTdenovo/$Organism/$Strain/nanopolish
mkdir -p $OutDir
# cat "" > $OutDir/"$Strain"_nanoplish.fa
InDir=$(dirname $Assembly)
NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_merge.py $InDir/*:*-*/*.fa > $OutDir/"$Strain"_nanoplish.fa

echo "" > tmp.txt
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $OutDir/"$Strain"_nanoplish.fa --out $OutDir/"$Strain"_nanoplish_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```

Quast and busco were run to assess the effects of nanopolish on assembly quality:

```bash
for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_albacore_v2/nanopolish/WT_albacore_v2_nanoplish_min_500bp_renamed.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
	# Quast
  OutDir=$(dirname $Assembly)
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
	# Busco
	BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
	OutDir=gene_pred/busco/$Organism/$Strain/assembly
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
	qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```



## Assemblies were polished using Pilon

 Although three libraries were available, the first contained a relatively small amount of data and was not used for correction.

<!-- As nanopolish was not run the assembly remained ralatively error prone. This
mena thtat the memorey overhead of pilon was very high. To reduce this, the
assembly was broken into contigs and an iteration of pilon run on each contig.
Following a single iteration of pilon correction this way, the assembly was merged once more and pilon run for a further x iterations on the combined assembly.

Split the assembly into contigs:
```bash
for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_albacore_v2/nanopolish/WT_albacore_v2_nanoplish_min_500bp_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)"/by_contig"
mkdir -p $OutDir
# Works only if fasta isnt wrapped:
cat $Assembly | split -l 2 - $OutDir/seq_
for File in $(ls $OutDir/* | grep -v 'fasta'); do
  mv $File $File.fasta
done
done
``` -->


```bash
for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_albacore_v2/nanopolish/WT_albacore_v2_nanoplish_min_500bp_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev | sed 's/_albacore_v2//g')
IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
echo $Strain
echo $Organism
#TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1 | tail -n1);
#TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1 | tail -n1);
TrimF2_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n2 | tail -n1);
TrimR2_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n2 | tail -n1);
TrimF3_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n3 | tail -n1);
TrimR3_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n3 | tail -n1);
#echo $TrimF1_Read
#echo $TrimR1_Read
echo $TrimF2_Read
echo $TrimR2_Read
echo $TrimF3_Read
echo $TrimR3_Read
OutDir=$(dirname $Assembly)
Iterations=10
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
#qsub $ProgDir/sub_pilon_3_libs.sh $Assembly $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir $Iterations
qsub $ProgDir/sub_pilon_2_libs.sh $Assembly $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir $Iterations
done
```


<!-- ```bash
for Assembly in $(ls assembly/canu-1.6/*/*/racon/*.fasta | grep 'WT' | grep 'round_10'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
echo $Strain
echo $Organism
#TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1 | tail -n1);
#TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1 | tail -n1);
TrimF2_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n2 | tail -n1);
TrimR2_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n2 | tail -n1);
TrimF3_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n3 | tail -n1);
TrimR3_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n3 | tail -n1);
#echo $TrimF1_Read
#echo $TrimR1_Read
echo $TrimF2_Read
echo $TrimR2_Read
echo $TrimF3_Read
echo $TrimR3_Read
OutDir=assembly/canu-1.6/$Organism/$Strain/polished_10
Iterations=10
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
#qsub $ProgDir/sub_pilon_3_libs.sh $Assembly $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir $Iterations
qsub $ProgDir/sub_pilon_2_libs.sh $Assembly $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir $Iterations
done
``` -->

Summarising assemblies using quast:
```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_albacore_v2/nanopolish/pilon_*.fasta | grep 'pilon_10'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```


checking using busco

```bash
for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_albacore_v2/nanopolish/pilon_*.fasta); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  # BuscoDB="Fungal"
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
  Prefix=$(echo $Assembly | rev | cut -d '/' -f1 | rev | sed 's/.fasta//g')
  OutDir=$(dirname $Assembly)
  OutDir=$OutDir/$Prefix
  echo "$Prefix"
  # OutDir=gene_pred/busco/$Organism/"$Strain"_pilon/assembly
  qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```


```bash
  for File in $(ls assembly/SMARTdenovo/F.venenatum/WT_albacore_v2/nanopolish/pilon_*/*/short_summary_*.txt); do  
  Strain=$(echo $File| rev | cut -d '/' -f5 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f6 | rev)
  Round=$(echo $File | rev | cut -d '/' -f3 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Round\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
```

```
F.venenatum	WT	pilon_1	3306	223	196	3725
F.venenatum	WT	pilon_2	3567	78	80	3725
F.venenatum	WT	pilon_3	3635	43	47	3725
F.venenatum	WT	pilon_4	3652	32	41	3725
F.venenatum	WT	pilon_5	3654	33	38	3725
F.venenatum	WT	pilon_6	3656	32	37	3725
F.venenatum	WT	pilon_7	3656	32	37	3725
F.venenatum	WT	pilon_8	3656	32	37	3725
F.venenatum	WT	pilon_9	3656	32	37	3725
F.venenatum	WT	pilon_10	3656	32	37	3725
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_albacore_v2/nanopolish/pilon_*.fasta | grep 'pilon_10'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
# OutDir=assembly/SMARTdenovo/$Organism/$Strain/nanopolish
OutDir=$(dirname $Assembly)
echo "" > tmp.txt
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/"$Strain"_pilon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```
<!--
### Hybrid assembly:

#### Hybrid assembly: Spades Assembly

```bash
  # for MinIonDat in $(ls raw_dna/minion/*/*/*_pass.fastq.gz); do
  # MinIonDat=$(ls raw_dna/minion/*/*/*_pass.fastq.gz)
  # mkdir tmp_assembly
  # cat $MinIonDat > tmp_assembly/miniondat.fastq.gz
  Reads1=$(ls qc_dna/minion/F.venenatum/WT/F.venenatum_WT_07-03-17_albacore_v2.02_trim.fastq.gz)
  Reads2=$(ls qc_dna/minion/F.venenatum/WT/F.venenatum_WT_18-07-17_albacore_v2.02_trim.fastq.gz)
  Organism=$(echo $Reads1 | head -n1 | rev | cut -f3 -d '/' | rev)
  Strain=$(echo $Reads1 | head -n1 | rev | cut -f2 -d '/' | rev)
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
  OutDir=assembly/spades_minion/$Organism/"$Strain"_albacore_v2
  # ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades/multiple_libraries
  # qsub $ProgDir/subSpades_3lib_minion.sh tmp_assembly/miniondat.fastq.gz $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_venenatum/assembly
  qsub $ProgDir/Fv_spades_hybrid.sh $Reads1 $Reads2 $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir
  # done
  rm -r tmp_assembly
```

```bash
  for Assembly in $(ls assembly/spades_minion/*/*/contigs.fasta |grep -v 'old' | grep '_albacore_v2'); do
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
  for Assembly in $(ls assembly/spades_minion/*/*/filtered_contigs/*_min_500bp.fasta |grep -v 'old' | grep '_albacore_v2'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    echo "$Organism - $Strain"
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
``` -->


# Merging Minion and Hybrid Assemblies
<!--
```bash
  for PacBioAssembly in $(ls assembly/canu-1.6/*/*/polished/*.fasta | grep -v -e '_nanopore' -e '_pass-fail'); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_*/$Organism/$Strain/contigs.fasta)
    OutDir=assembly/merged_canu_spades/$Organism/$Strain
    AnchorLength=500000
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
  done
``` -->
<!--
```bash
  for PacBioAssembly in $(ls assembly/canu-1.6/*/*/polished/*.fasta | grep -v -e '_nanopore' -e '_pass-fail'); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_*/$Organism/$Strain/contigs.fasta)
    OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_100k_anchor
    AnchorLength=100000
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
  done
``` -->
<!--
```bash
  for PacBioAssembly in $(ls assembly/canu-1.6/F.venenatum/WT/polished_10/pilon_*.fasta | grep 'pilon_10'); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_*/$Organism/$Strain/contigs.fasta)
    # OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_spades_first_corrected
    OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_spades_first_100k
    AnchorLength=100000
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $HybridAssembly $PacBioAssembly $OutDir $AnchorLength
  done
```

```bash
  for PacBioAssembly in $(ls assembly/canu-1.6/F.venenatum/WT/polished_10/pilon_*.fasta | grep 'pilon_10'); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_*/$Organism/$Strain/contigs.fasta)
    OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_minion_first_100k
    # OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_minion_first_corrected
    AnchorLength=100000
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
  done
```

```bash
  for PacBioAssembly in $(ls assembly/canu-1.4/F.venenatum/WT/polished_10/pilon_*.fasta | grep 'pilon_10'); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_*/$Organism/$Strain/contigs.fasta)
    AnchorLength=100000
    OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_minion_first_100k_canu_1.4
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
    OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_spades_first_100k_canu_1.4
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $HybridAssembly $PacBioAssembly $OutDir $AnchorLength
  done
```

Checking assembly quality

```bash
# for Assembly in $(ls assembly/merged_canu_spades/*/*/merged.fasta | grep -e '_spades_first_corrected' -e 'hybrid_first_corrected' | grep -v 'old' | grep 'hybrid_first_corrected'); do
for Assembly in $(ls assembly/merged_canu_spades/*/*_100k_canu_1.4/merged.fasta); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

checking using busco

```bash
#for Assembly in  $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa); do
for Assembly in $(ls assembly/merged_canu_spades/*/*/merged.fasta | grep -e '_spades_first_corrected' -e 'hybrid_first_corrected' | grep -v 'old' | grep 'hybrid_first_corrected'); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  # BuscoDB="Fungal"
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
  OutDir=gene_pred/busco/$Organism/$Strain/assembly
  qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
  for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt); do  
    echo $File;
    cat $File | grep -e '(C)' -e '(F)' -e '(M)' -e 'Total';
  done
```
-->

<!--
This merged assembly was polished using Pilon

```bash
for Assembly in $(ls assembly/merged_canu_spades/*/*/merged.fasta | grep '_spades_first_corrected' | grep -v 'old'); do
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev | cut -f1 -d '_')
IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
echo $Strain
echo $Organism
#TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1 | tail -n1);
#TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1 | tail -n1);
TrimF2_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n2 | tail -n1);
TrimR2_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n2 | tail -n1);
TrimF3_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n3 | tail -n1);
TrimR3_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n3 | tail -n1);
#echo $TrimF1_Read
#echo $TrimR1_Read
echo $TrimF2_Read
echo $TrimR2_Read
echo $TrimF3_Read
echo $TrimR3_Read
Iterations=3
OutDir=$(dirname $Assembly)"/polished$Iterations"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
qsub $ProgDir/sub_pilon_2_libs.sh $Assembly $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir $Iterations
done
```

## Editing contigs

Contigs were renamed in accordance with ncbi recomendations

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  for Assembly in $(ls assembly/merged_canu_spades/F.venenatum/WT_spades_first_corrected/polished*/pilon_3.fasta); do
    TrimInstructions=tmp.txt
    printf "" > $TrimInstructions
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev | cut -f1 -d '_')
    OutDir=$(dirname $Assembly)
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_merged_polished.fasta --coord_file $TrimInstructions
  done
```


Checking assembly quality

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/merged_canu_spades/F.venenatum/WT_spades_first_corrected/polished*/*_merged_polished.fasta); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

checking using busco

```bash
for Assembly in $(ls assembly/merged_canu_spades/F.venenatum/WT_spades_first_corrected/polished*/*_merged_polished.fasta); do
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev | cut -f1 -d '_')
  echo "$Organism - $Strain"
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  # BuscoDB="Fungal"
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
  OutDir=gene_pred/busco/$Organism/$Strain/assembly
  qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
  for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt); do  
    echo $File;
    cat $File | grep -e '(C)' -e '(F)' -e '(M)' -e 'Total';
  done
```
-->

# Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_albacore_v2/nanopolish/WT_albacore_v2_pilon_min_500bp_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev | cut -f1 -d '_')
echo "$Organism - $Strain"
OutDir=repeat_masked/$Organism/"$Strain"_minion/minion_submission
qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
done
```


The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and hardmasked files.


```bash
for File in $(ls repeat_masked/F.venenatum/WT_minion/minion_submission/*_contigs_softmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/F.venenatum/WT_minion/minion_submission/*_contigs_softmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```


```bash
for RepDir in $(ls -d repeat_masked/F.venenatum/WT_minion/minion_submission); do
Strain=$(echo $RepDir | rev | cut -f2 -d '/' | rev)
Organism=$(echo $RepDir | rev | cut -f3 -d '/' | rev)  
RepMaskGff=$(ls $RepDir/*_contigs_hardmasked.gff)
TransPSIGff=$(ls $RepDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
# printf "The number of bases masked by RepeatMasker:\t"
RepMaskerBp=$(sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
# printf "The number of bases masked by TransposonPSI:\t"
TpsiBp=$(sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
# printf "The total number of masked bases are:\t"
Total=$(cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
printf "$Organism\t$Strain\t$RepMaskerBp\t$TpsiBp\t$Total\n"
done
```

```bash
  F.venenatum	WT_minion	476671	128276	602385
```


# Gene Prediction

Gene prediction followed steps:
Gene model training
		- Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline
	Gene prediction
		- Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.


#### Aligning


```bash
  for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/minion_submission/*_contigs_unmasked.fa); do
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
OutDir=alignment/star/F.venenatum/WT_minion/treatment/concatenated
mkdir -p $OutDir
# for File in $(ls alignment/star/F.venenatum/WT_minion/treatment/*/*.sortedByCoord.out.bam); do
# Prefix=$(basename $File | sed 's/.bam//g')
# echo $Prefix
# samtools sort -o $File $Prefix > $OutDir/"$Prefix"_sorted.bam
# done
BamFiles=$(ls alignment/star/F.venenatum/WT_minion/treatment/*/*.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
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


Results of web-annotation of gene clusters within the assembly were downloaded to
the following directories:

```bash
  for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep 'WT_ncbi'); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=gene_pred/secondary_metabolites/antismash/$Organism/$Strain
    mkdir -p $OutDir
  done
```

```bash
  for Zip in $(ls gene_pred/secondary_metabolites/antismash/*/*/*.zip); do
    OutDir=$(dirname $Zip)
    unzip -d $OutDir $Zip
  done
```


```bash
  for AntiSmash in $(ls gene_pred/secondary_metabolites/antismash/*/*/*/*.final.gbk); do
    Organism=$(echo $AntiSmash | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $AntiSmash | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/secondary_metabolites/antismash/$Organism/$Strain
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
    $ProgDir/antismash2gff.py --inp_antismash $AntiSmash > $OutDir/"$Strain"_secondary_metabolite_regions.gff
    printf "Number of clusters detected:\t"
    cat $OutDir/"$Strain"_secondary_metabolite_regions.gff | grep 'antismash_cluster' | wc -l
    # GeneGff=gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3
    GeneGff=gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3
    bedtools intersect -u -a $GeneGff -b $OutDir/"$Strain"_secondary_metabolite_regions.gff > $OutDir/metabolite_cluster_genes.gff
    cat $OutDir/metabolite_cluster_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > $OutDir/metabolite_cluster_gene_headers.txt
    printf "Number of predicted proteins in clusters:\t"
    cat $OutDir/metabolite_cluster_gene_headers.txt | wc -l
    printf "Number of predicted genes in clusters:\t"
    cat $OutDir/metabolite_cluster_genes.gff | grep -w 'gene' | wc -l
  done
```
<!--
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
  for Gff in $(ls gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.gff3); do
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
