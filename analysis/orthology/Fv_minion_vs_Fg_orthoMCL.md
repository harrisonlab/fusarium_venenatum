# Fv vs Fg

Orthology analysis between Fv and Fg

Using most up-to-date Fg gene models:

```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium_venenatum
  cd $ProjDir
  IsolateAbrv=Fv_minion_vs_Fg_JGI
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
```

## 2.1 Format fasta files


### for Fv WT (strain name A3/5)
```bash
  Taxon_code=A3_5
  Fasta_file=$(ls assembly/external_group/F.graminearum/PH1/pep/Fusarium_graminearum.RR1.pep.all.fa)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for Fg PH-1

```bash
  Taxon_code=FusGr1
  Fasta_file=$(ls assembly/external_group/F.graminearum/Fusgr1/Fusgr1_GeneCatalog_proteins_20110524.aa.fasta)
  Id_field=4
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```



## 2.2 Filter proteins into good and poor sets.

```bash
  Input_dir=$WorkDir/formatted
  Min_length=10
  Max_percent_stops=20
  Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
  Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
  orthomclFilterFasta $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
```
<!--
## 2.3.a Perform an all-vs-all blast of the proteins

```bash
  BlastDB=$WorkDir/blastall/$IsolateAbrv.db

  makeblastdb -in $Good_proteins_file -dbtype prot -out $BlastDB
  BlastOut=$WorkDir/all-vs-all_results.tsv
  mkdir -p $WorkDir/splitfiles

  SplitDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
  $SplitDir/splitfile_500.py --inp_fasta $Good_proteins_file --out_dir $WorkDir/splitfiles --out_base goodProteins

  ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology  
  for File in $(find $WorkDir/splitfiles); do
    Jobs=$(qstat | grep 'blast_500' | grep 'qw' | wc -l)
    while [ $Jobs -gt 1 ]; do
      sleep 3
      printf "."
      Jobs=$(qstat | grep 'blast_500' | grep 'qw' | wc -l)
    done
    printf "\n"
    echo $File
    BlastOut=$(echo $File | sed 's/.fa/.tab/g')
    qsub $ProgDir/blast_500.sh $BlastDB $File $BlastOut
  done
```

## 2.3.b Merge the all-vs-all blast results  
```bash  
  MergeHits="$IsolateAbrv"_blast.tab
  printf "" > $MergeHits
  for Num in $(ls $WorkDir/splitfiles/*.tab | rev | cut -f1 -d '_' | rev | sort -n); do
    File=$(ls $WorkDir/splitfiles/*_$Num)
    cat $File
  done > $MergeHits
```
-->



## using orthofinder

```bash
qlogin -pe smp 8

#16 threads used
ProjDir=/home/groups/harrisonlab/project_files/fusarium_venenatum
cd $ProjDir
IsolateAbrv=Fv_minion_vs_Fg_JGI
WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
orthofinder -f $WorkDir/formatted -t 3 -a 3
```

orthofinder results:

```
OrthoFinder assigned 157358 genes (99.2% of total) to 14200 orthogroups. Fifty percent of all genes were in orthogroups
with 12 or more genes (G50 was 12) and were contained in the largest 6182 orthogroups (O50 was 6182). There were 10919
orthogroups with all species present and 10285 of these consisted entirely of single-copy genes.
```


## 2.4 Perform ortholog identification

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  MergeHits="$IsolateAbrv"_blast.tab
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  qsub $ProgDir/qsub_orthomcl.sh $MergeHits $GoodProts 5
```

## 2.5.b Plot venn diagrams:

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/venn_diagrams
  $ProgDir/venn_diag_2_way.r --inp $WorkDir/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/"$IsolateAbrv"_orthogroups.pdf
```

Output was a pdf file of the venn diagram.

The following additional information was also provided. The format of the
following lines is as follows:

Isolate name (total number of orthogroups)
number of unique singleton genes
number of unique groups of inparalogs


```
[1] "A3_5 (10933)"
[1] 1490
[1] 59
[1] "FusGr1 (11300)"
[1] 1892
[1] 24
```


#### 2.6) Extracting fasta files for all orthogroups

```bash
WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
GoodProt=$WorkDir/goodProteins/goodProteins.fasta
OutDir=$WorkDir/orthogroups_fasta
mkdir -p $OutDir
$ProgDir/orthoMCLgroups2fasta.py --orthogroups $WorkDir/"$IsolateAbrv"_orthogroups.txt --fasta $GoodProt --out_dir $OutDir > $OutDir/extractionlog.txt
```
