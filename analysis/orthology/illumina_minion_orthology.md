# Fv vs Fg

Orthology analysis between Fv illimina and Fv minion


# 1. Using old Fg gene models:

```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium_venenatum
  cd $ProjDir
  IsolateAbrv=MiSeq_vs_MinION
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
```

## 1.1 Format fasta files


### for Fv WT illumina MiSeq (strain name A3/5)
```bash
  Taxon_code=MiSeq
  Fasta_file=$(ls gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for Fv WT ONT minion (strain name A3/5)

```bash
Taxon_code=MinION
Fasta_file=$(ls gene_pred/final/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta)
Id_field=1
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```



## 1.2 Filter proteins into good and poor sets.

```bash
  Input_dir=$WorkDir/formatted
  Min_length=10
  Max_percent_stops=20
  Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
  Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
  orthomclFilterFasta $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
```

## 1.3.a Perform an all-vs-all blast of the proteins

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

## 1.3.b Merge the all-vs-all blast results  
```bash  
  MergeHits="$IsolateAbrv"_blast.tab
  printf "" > $MergeHits
  for Num in $(ls $WorkDir/splitfiles/*.tab | rev | cut -f1 -d '_' | rev | sort -n); do
    File=$(ls $WorkDir/splitfiles/*_$Num)
    cat $File
  done > $MergeHits
```

## 1.4 Perform ortholog identification

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  MergeHits="$IsolateAbrv"_blast.tab
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  qsub $ProgDir/qsub_orthomcl.sh $MergeHits $GoodProts 5
```


Also try using orthofinder

```bash
  qlogin -pe smp 8

  #16 threads used
  ProjDir=/home/groups/harrisonlab/project_files/fusarium_venenatum
  cd $ProjDir
  IsolateAbrv=MiSeq_vs_MinION
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  orthofinder -f $WorkDir/formatted -t 8 -a 8
```

orthofinder results:

```
OrthoFinder assigned 157080 genes (99.2% of total) to 14190 orthogroups. Fifty percent of all genes were in orthogroups
with 12 or more genes (G50 was 12) and were contained in the largest 6161 orthogroups (O50 was 6161). There were 10661
orthogroups with all species present and 10011 of these consisted entirely of single-copy genes.
```

output files are in:
```bash
ls $WorkDir/formatted/Results_Apr10
```

## 1.5.b Plot venn diagrams:

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
[1] "A3_5 (10897)"
[1] 1323
[1] 57
[1] "PH1 (12126)"
[1] 2581
[1] 28
```


#### 1.6) Extracting fasta files for all orthogroups

```bash
WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
GoodProt=$WorkDir/goodProteins/goodProteins.fasta
OutDir=$WorkDir/orthogroups_fasta
mkdir -p $OutDir
$ProgDir/orthoMCLgroups2fasta.py --orthogroups $WorkDir/"$IsolateAbrv"_orthogroups.txt --fasta $GoodProt --out_dir $OutDir > $OutDir/extractionlog.txt
```



# 2. Using most up-to-date Fg gene models:

```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium_venenatum
  cd $ProjDir
  IsolateAbrv=Fv_vs_Fg_JGI
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
  Fasta_file=$(ls gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.pep.fasta)
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

## 2.4 Perform ortholog identification

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  MergeHits="$IsolateAbrv"_blast.tab
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  qsub $ProgDir/qsub_orthomcl.sh $MergeHits $GoodProts 5
```
<!--
## 3.5.a Manual identification of numbers of orthologous and unique genes


```bash
  for num in 1; do
    echo "The number of ortholog groups found in pathogen but absent in non-pathogens is:"
    cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'PG|' -e 'A13|' -e 'fo47|' -e 'CB3|' | grep 'Fus2|' | grep '125|' | grep 'A23|' |  wc -l
    echo "The number of ortholog groups unique to pathogens are:"
    cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'PG|' -e 'A13|' -e 'fo47|' -e 'CB3|' | grep 'Fus2|' | grep '125|' | grep 'A23|' | wc -l
    echo "The number of ortholog groups unique to non-pathogens are:"
    cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'Fus2|' -e '125|' -e 'A23|' | grep 'A13|' | grep 'A28|' | grep 'PG|' | grep 'fo47|' | grep 'CB3|' | wc -l
    echo "The number of ortholog groups common to all F. oxysporum isolates are:"
    cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep 'A28|' | grep 'PG|' | grep 'A13|' | grep 'fo47|' | grep 'CB3|' | grep '4287|' |wc -l
  done
```

```
  The number of ortholog groups found in pathogen but absent in non-pathogens is:
  277
  The number of ortholog groups unique to pathogens are:
  277
  The number of ortholog groups unique to non-pathogens are:
  64
  The number of ortholog groups common to all F. oxysporum isolates are:
  10396
```

The number of ortholog groups shared between FoN and FoL was identified:

```bash
  echo "The number of ortholog groups common to FoN and FoL are:"
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'N139|' | grep '4287|' | wc -l
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'PG|' -e 'A13|' -e 'fo47|' -e 'CB3' | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep '4287|' |  wc -l
```

```
  The number of ortholog groups common to FoC and FoL are:
  10629
  37
``` -->

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
