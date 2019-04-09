# Promoter analysis

Commands to identify regulatory elements in promoters of gene clusters

## Extract promotor regions:

Promotor regions upstream of all genes were extracted:

```bash

for GeneGff in $(ls gene_pred/final/F.*/*/*/final_genes_appended_renamed.gff3 | grep -w 'WT'); do
Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep 'ncbi')
OutDir=analysis/meme/promotor_regions/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/RNAseq/P414/motif_discovery
$ProgDir/extract_promotor.pl --fasta $Assembly --gff $GeneGff --prefix $OutDir/${Strain} --ranges 1:100 101:200 201:300 301:400 401:500
done

ls $OutDir
```

# Tri5 analysis:

Two domains are already characteristed in the Tri cluster. Analysis of this
region allows these to be used as positive controls for the analysis.

Motifs:
* YNAGGCC
* GTGA/TCAC(8bp then repeat) Zn-finger binding


## Extract cluster intergenic regions

The Tri cluster has been identified as secondary metabolite cluster 15, with the
 known cluster genes identified as:
```
Fgtyr - g3123
Tri8 - g3124
(Missing) - g3125
Tri3 - g3126
Tri4 - g3127
Tri6 - g3128
Tri5 - g3129
Tri10 - g3130
Tri9 - (not predicted)
Tri11 - g3131
Tri12 - g3132
Tri13 - g3133
Tri14 - g3134
Fgest - g3135
```

These Tri5 genes were extracted:

```bash

OutDir=analysis/meme/promotor_regions/F.venenatum/WT/tri
mkdir -p $OutDir

TriGenes="g3123.t1 g3124.t1 g3125.t1 g3126.t1 g3127.t1 g3128.t1 g3129.t1 g3130.t1 g3131.t1 g3132.t1 g3133.t1 g3134.t1 g3135.t1"
printf "$TriGenes" | sed "s/ /\n/g" | sed "s/.t.//g" > $OutDir/tri_headers.txt

# Create fasta files of RxLR upstream regions
for Upstream in $(ls analysis/meme/promotor_regions/F.venenatum/WT/*.upstream*.fasta); do
  Region=$(basename ${Upstream%.fasta} | sed 's/promotor_regions.upstream//g')
  mkdir $OutDir/$Region
  RegionPromotors=$OutDir/$Region/F.venenatum_WT_${Region}_promotors.fa
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $Upstream --headers $OutDir/tri_headers.txt > $RegionPromotors
done
```


```bash
# Run MEME on  promotor regions

screen -a

qlogin -pe smp 4
cd /home/groups/harrisonlab/project_files/fusarium_venenatum
for Query in $(ls analysis/meme/promotor_regions/F.venenatum/WT/tri/_promotors.fa); do
  Region=$(basename ${Query%.fa} | sed 's/promotor_regions.upstream//g')
  echo $Region
  OutDir=$(dirname $Query)
  mkdir -p $OutDir/meme
  meme $Query -p 4 -dna -oc $OutDir/meme -nostatus -mod zoops -nmotifs 5 -minw 6 -maxw 20 -objfun classic -revcomp -markov_order 0
  ls $OutDir/meme/meme.txt
  mast $OutDir/meme/meme.xml $Query -oc $OutDir/meme -nostatus
  mv $OutDir/meme/mast.txt $OutDir/meme/${Region}_mast.txt
  mv $OutDir/meme/mast.html $OutDir/meme/${Region}_mast.html
done

# Also initial run of DREME on the webserver:
# dreme -verbosity 1 -oc . -dna -p P414_RxLR_promotors.fa -t 18000 -e 0.05 -dfile description
# qlogin -pe smp 4
# cd /data/scratch/armita/idris
# OutDir=analysis/meme/RxLR
# RxlrPromotors=$(ls $OutDir/P414_RxLR_promotors.fa)
# mkdir -p $OutDir/dreme
# dreme -verbosity 1 -oc $OutDir/dreme -dna -p $RxlrPromotors -t 18000 -e 0.05

for Query in $(ls analysis/meme/RxLR/*/P414_RxLR_*_promotors.fa); do
  Region=$(basename ${Query%.fa} | sed 's/promotor_regions.upstream//g')
  echo $Region
  OutDir=$(dirname $Query)
  mkdir -p $OutDir/dreme
  dreme -verbosity 1 -oc $OutDir/dreme -dna -p $Query -e 0.05
  ls $OutDir/dreme/dreme.txt
  mast $OutDir/dreme/dreme.xml $Query -oc $OutDir/dreme -nostatus
  mv $OutDir/dreme/mast.txt $OutDir/dreme/${Region}_mast.txt
  mv $OutDir/dreme/mast.html $OutDir/dreme/${Region}_mast.html
done
```





Cluster15 genes were extracted


Blast searching identified Fv g3129 as the homolog of Fv tri5

It didnt show expression on myro media in the RNAseq analysis
```bash
cat /home/deakig/projects/quorn/DGE/all.merged.tsv | grep -e 'rowname' -e 'g3129' |cut -f9
```
```
FC_MWT
-2.43870548122197
```


```bash


```
