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
for Query in $(ls analysis/meme/promotor_regions/F.venenatum/WT/tri/*/*_promotors.fa); do
  Region=$(basename ${Query%.fa} | sed 's/.*upstream//g' | sed "s/_prom.*//g")
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

for Query in $(ls analysis/meme/promotor_regions/F.venenatum/WT/tri/*/*_promotors.fa); do
  Region=$(basename ${Query%.fa} | sed 's/.*upstream//g' | sed "s/_prom.*//g")
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

## 1:1000 bp analysis:

### Extract promotor regions:

Promotor regions upstream of all genes were extracted:

```bash
for GeneGff in $(ls gene_pred/final/F.*/*/*/final_genes_appended_renamed.gff3 | grep -w 'WT'); do
Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep 'ncbi')
OutDir=analysis/meme/promotor_regions/1-1000/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/RNAseq/P414/motif_discovery
$ProgDir/extract_promotor.pl --fasta $Assembly --gff $GeneGff --prefix $OutDir/${Strain} --ranges 1:1000
done

ls $OutDir
```

Tri5 genes were extracted:

```bash
OutDir=analysis/meme/promotor_regions/1-1000/F.venenatum/WT/tri
mkdir -p $OutDir

TriGenes="g3123.t1 g3124.t1 g3125.t1 g3126.t1 g3127.t1 g3128.t1 g3129.t1 g3130.t1 g3131.t1 g3132.t1 g3133.t1 g3134.t1 g3135.t1"
printf "$TriGenes" | sed "s/ /\n/g" | sed "s/.t.//g" > $OutDir/tri_headers.txt

# Create fasta files of RxLR upstream regions
for Upstream in $(ls analysis/meme/promotor_regions/1-1000/F.venenatum/WT/*.upstream*.fasta); do
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
for Query in $(ls analysis/meme/promotor_regions/1-1000/F.venenatum/WT/tri/*/*_promotors.fa); do
  Region=$(basename ${Query%.fa} | sed 's/.*upstream//g' | sed "s/_prom.*//g")
  echo $Region
  OutDir=$(dirname $Query)
  mkdir -p $OutDir/meme
  meme $Query -p 4 -dna -oc $OutDir/meme -nostatus -mod zoops -nmotifs 5 -minw 6 -maxw 20 -objfun classic -revcomp -markov_order 0
  ls $OutDir/meme/meme.txt
  mast $OutDir/meme/meme.xml $Query -oc $OutDir/meme -nostatus
  mv $OutDir/meme/mast.txt $OutDir/meme/${Region}_mast.txt
  mv $OutDir/meme/mast.html $OutDir/meme/${Region}_mast.html
done

for Query in $(ls analysis/meme/promotor_regions/1-1000/F.venenatum/WT/tri/*/*_promotors.fa); do
  Region=$(basename ${Query%.fa} | sed 's/.*upstream//g' | sed "s/_prom.*//g")
  echo $Region
  OutDir=$(dirname $Query)
  mkdir -p $OutDir/dreme
  dreme -verbosity 1 -oc $OutDir/dreme -dna -p $Query -e 0.05
  ls $OutDir/dreme/dreme.txt
  mast $OutDir/dreme/dreme.xml $Query -oc $OutDir/dreme -nostatus
  mv $OutDir/dreme/mast.txt $OutDir/dreme/${Region}_mast.txt
  mv $OutDir/dreme/mast.html $OutDir/dreme/${Region}_mast.html
done


for Query in $(ls analysis/meme/promotor_regions/1-1000/F.venenatum/WT/tri/*/*_promotors.fa); do
  Region=$(basename ${Query%.fa} | sed 's/.*upstream//g' | sed "s/_prom.*//g")
  AllGenes=$(ls analysis/meme/promotor_regions/1-1000/F.venenatum/WT/*.upstream*.fasta)
  NumQueries=$(cat $Query | grep '>' | wc -l)
  RANDOM=$(date +%N | sed s/...$//)
  Control=${Query%.*})_random_$NumQueries.fasta
  /home/sobczm/bin/meme_4.11.2/bin/fasta-subsample -seed $RANDOM $AllGenes $NumQueries > $Control
  echo $Region
  OutDir=$(dirname $Query)
  mkdir -p $OutDir/dreme
  dreme -verbosity 1 -oc $OutDir/dreme -n $Control -dna -p $Query -e 0.05
  ls $OutDir/dreme/dreme.txt
  mast $OutDir/dreme/dreme.xml $Query -oc $OutDir/dreme -nostatus
  mv $OutDir/dreme/mast.txt $OutDir/dreme/${Region}_mast.txt
  mv $OutDir/dreme/mast.html $OutDir/dreme/${Region}_mast.html
done
```


## Running Maria's pipeline

https://github.com/harrisonlab/popgen/tree/master/clock

```
Script execution order:

clock_ortho.sh && clock_ortho_andy.sh
extract_promoter.sh
clock_motif_discovery.sh
clock_motif_discovery_cont.sh
clock_motif_discovery2.sh
clock_dn_ds_analysis.sh
dn_ds_analysis_pairwise_verticillium.sh
domain_analysis.sh
```

### Enrichment of known motifs

using meme ame

```bash
  UpstreamFa=$(ls gene_pred/final/F.venenatum/WT_ncbi/final/final_genes_combined.upstream3000.fasta)
  OutDir=analysis/meme/promotor_regions/F.venenatum/WT/tri
  mkdir -p $OutDir
  cat $UpstreamFa | sed -n "/g3123/,/g3136/p;/g3136/q" | head -n-1 > $OutDir/tri_cluster_upstream3000.fasta

  #Create a random sample of a 100 promoter control sequences

  scripts=/home/sobczm/bin/popgen/clock/motif_discovery
  # qsub $scripts/sub_fasta_subsample.sh gene_pred/final/F.venenatum/WT_ncbi/final/final_genes_combined.upstream3000.fasta 100
  qsub $scripts/sub_fasta_subsample.sh gene_pred/final/F.venenatum/WT_ncbi/final/final_genes_combined.upstream3000.fasta 13

  OutDir=analysis/meme/promotor_regions/F.venenatum/WT/tri
  # qsub $scripts/sub_ame.sh $OutDir/tri_cluster_upstream3000.fasta final_genes_combined.upstream3000_random_100.fasta min_ace YNAGGCC
  qsub $scripts/sub_ame.sh $OutDir/tri_cluster_upstream3000.fasta final_genes_combined.upstream3000_random_13.fasta Zn_finger GTGA
  mv *_vs_Zn_finger $OutDir/.
```


motif searches using motif scanning:

```bash
#As no motif enrichment can be conducted that way, motif scanning will be carried out.
# cd $dna/promoters/extended/ace
scripts=/home/sobczm/bin/popgen/clock/motif_discovery
# qsub $scripts/sub_fimo.sh $OutDir/tri_cluster_upstream3000.fasta Zn_finger GTGA
# mv tri_cluster_upstream3000_vs_Zn_finger $OutDir/.
qsub $scripts/sub_fimo.sh $OutDir/tri_cluster_upstream3000.fasta Tri6_binding YNAGGCC
mv tri_cluster_upstream3000_vs_Tri6_binding $OutDir/.
cat $OutDir/tri_cluster_upstream3000_vs_Tri6_binding/fimo.txt | tail -n+2 | wc -l
cat $OutDir/tri_cluster_upstream3000_vs_Tri6_binding/fimo.txt | tail -n+2 | cut -f2 | sort | uniq -c

scripts=/home/sobczm/bin/popgen/clock/motif_discovery
qsub $scripts/sub_fasta_subsample.sh gene_pred/final/F.venenatum/WT_ncbi/final/final_genes_combined.upstream3000.fasta 13
mv final_genes_combined.upstream3000_random_13.fasta $OutDir/.
qsub $scripts/sub_fimo.sh $OutDir/final_genes_combined.upstream3000_random_13.fasta Tri6_binding YNAGGCC
mv final_genes_combined.upstream3000_random_13_vs_Tri6_binding $OutDir/.
cat $OutDir/final_genes_combined.upstream3000_random_13_vs_Tri6_binding/fimo.txt | tail -n+2 | wc -l
cat $OutDir/final_genes_combined.upstream3000_random_13_vs_Tri6_binding/fimo.txt | tail -n+2 | cut -f2 | sort | uniq -c
```
```
4103

273 g3123_upstream3000
401 g3124_upstream3000
398 g3125_upstream3000
342 g3126_upstream3000
354 g3127_upstream3000
297 g3128_upstream3000
283 g3129_upstream3000
322 g3130_upstream3000
266 g3131_upstream3000
341 g3132_upstream3000
258 g3133_upstream3000
280 g3134_upstream3000
288 g3135_upstream3000

3930

325 g10111_upstream3000
327 g10414_upstream3000
278 g10565_upstream3000
303 g12182_upstream3000
294 g14911_upstream3000
346 g16730_upstream3000
283 g18968_upstream3000
272 g2035_upstream3000
328 g2249_upstream3000
291 g5643_upstream3000
293 g6086_upstream3000
314 g6708_upstream3000
276 g9898_upstream3000
```

### DREME and GLAM2
```bash
  TriPromoters=$(ls $OutDir/tri_cluster_upstream3000.fasta)
  Control=$(ls $OutDir/final_genes_combined.upstream3000_random_13.fasta)
  scripts=/home/sobczm/bin/popgen/clock/motif_discovery
  qsub $scripts/sub_dreme.sh $TriPromoters $Control
  qsub $scripts/sub_glam.sh $TriPromoters

  mv tri_cluster_upstream3000_dreme $OutDir/.
  mv tri_cluster_upstream3000_glam $OutDir/.
```



<!-- ACE element in ccg2

The core ACE sequence binding site (Bell-Pedersen 2001), from -107,
similar a bit to a sequence in other ccgs (Corrego 2003)
AACTTGGCCAAGTT

Use FIMO

```bash
qsub $scripts/sub_fimo.sh Neurospora_crassa.NC12.dna_rm.toplevel_promoters_1000.fasta \
core_ace AACTTGGCCAAGTT
``` -->

```bash
###############Motif enrichment testing
# A) ACE motif containing (TCTTGGCA)
# B) Clockbox motif containing (CGAT(N)CCGCT)

input=/home/sobczm/popgen/clock/DNA_genomes/promoters
scripts=/home/sobczm/bin/popgen/clock/motif_discovery
meme=/home/sobczm/bin/meme_4.11.2/bin

#Create a random sample of a 100 promoter control sequences
qsub $scripts/sub_fasta_subsample.sh $dna/Neurospora_crassa.NC12.dna_rm.toplevel_promoters_1000.fasta 37
qsub $scripts/sub_fasta_subsample.sh $dna/Neurospora_crassa.NC12.dna_rm.toplevel_promoters_2000.fasta 100

#Test against the enrichment of motifs
#!!Warning!! it does not allow for motifs of variable length, like the cbox
#Create a random sample of a 100 promoter control sequences
qsub $scripts/sub_fasta_subsample.sh $dna/Neurospora_crassa.NC12.dna_rm.toplevel_promoters_1000.fasta 37
cd $dna/promoters/extended/ace
qsub $scripts/sub_ame.sh ace_1000_nc.fasta Neurospora_crassa.NC12.dna_rm.toplevel_promoters_1000_random_37.fasta \
min_ace TCTTGGCA
qsub $scripts/sub_fasta_subsample.sh $dna/Neurospora_crassa.NC12.dna_rm.toplevel_promoters_2000.fasta 37
qsub $scripts/sub_ame.sh ace_2000_nc.fasta Neurospora_crassa.NC12.dna_rm.toplevel_promoters_2000_random_37.fasta \
min_ace TCTTGGCA

cd $dna/promoters/extended/cbox
qsub $scripts/sub_ame.sh cbox_1000_nc.fasta Neurospora_crassa.NC12.dna_rm.toplevel_promoters_1000_random.fasta \
min_cbox1 CGAT
qsub $scripts/sub_ame.sh cbox_1000_nc.fasta Neurospora_crassa.NC12.dna_rm.toplevel_promoters_1000_random.fasta \
min_cbox2 CCGCT
#As no motif enrichment can be conducted that way, motif scanning will be carried out.
cd $dna/promoters/extended/ace
qsub $scripts/sub_fimo.sh ace_2000_nc.fasta min_ace TCTTGGCA

#Found some motifs on both strands but mostly imperfect, and even inside the exons!
#Cannot replicate the results from Correa et al. (2003) even when checking for the sequences manually.
```


# Cluster 15

## Cluster15 genes were extracted


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
