# Fv vs Fg

Orthology analysis between Fv and Fg

Using most up-to-date Fg gene models:

```bash
ProjDir=/projects/fusarium_venenatum
cd $ProjDir
IsolateAbrv=Fv_vs_Fg_JGI
WorkDir=analysis/orthology/$IsolateAbrv
mkdir -p $WorkDir
mkdir -p $WorkDir/formatted
mkdir -p $WorkDir/goodProteins
mkdir -p $WorkDir/badProteins  
```

## 2.1 Format fasta files


### for Fv WT minion genome (strain name A3/5)
```bash
Taxon_code=WT_M
Fasta_file=$(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta)
Id_field=1
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for Fv WT MiSeq genome (strain name A3/5)
```bash
  Taxon_code=WT
  Fasta_file=$(ls gene_pred/codingquarry/F.venenatum/WT/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for Fg PH-1

```bash
  Taxon_code=FusGr1
  Fasta_file=$(ls ../oldhome/groups/harrisonlab/project_files/fusarium_venenatum/assembly/external_group/F.graminearum/Fusgr1/Fusgr1_GeneCatalog_proteins_20110524.aa.fasta)
  Id_field=4
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

## using orthofinder

```bash

srun --partition long --mem 200G --cpus-per-task 20 --pty bash

#16 threads used

IsolateAbrv=Fv_vs_Fg_JGI
WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
orthofinder -f formatted -t 3 -a 6

IsolateAbrv=Fv_MiSeq_vs_Fg_JGI
WorkDir=analysis/orthology/$IsolateAbrv
mkdir -p $WorkDir
mkdir -p $WorkDir/formatted
mkdir -p $WorkDir/goodProteins
mkdir -p $WorkDir/badProteins  

### for Fv WT MiSeq genome (strain name A3/5)
```bash
Taxon_code=WT
Fasta_file=$(ls gene_pred/codingquarry/F.venenatum/WT/final/final_genes_appended_renamed.pep.fasta)
Id_field=1
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for Fg PH-1

```bash
Taxon_code=FusGr1
Fasta_file=$(ls ../oldhome/groups/harrisonlab/project_files/fusarium_venenatum/assembly/external_group/F.graminearum/Fusgr1/Fusgr1_GeneCatalog_proteins_20110524.aa.fasta)
Id_field=4
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```



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