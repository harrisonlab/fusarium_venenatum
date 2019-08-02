

## Stand alone run of CASSIS


```bash
screen -a
ssh compute02

conda active meme-v4

WorkDir=~/tmp/cassis_Fv
mkdir $WorkDir
cd $WorkDir
OldProjDir=/oldhpc/home/groups/harrisonlab/project_files/fusarium_venenatum
Assembly=$(ls $OldProjDir/repeat_masked/F.venenatum/WT/illumina_assembly_ncbi/WT_contigs_softmasked_repeatmasker_TPSI_appended.fa)
Genes=$(ls $OldProjDir/gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.gff3)
# Genes=$(ls $OldProjDir/gene_pred/braker/F.venenatum/WT_UTR/*/augustus_extracted.gff)
Interpro=$(ls $OldProjDir/gene_pred/interproscan/F.venenatum/WT/WT_interproscan.tsv)

cp $Interpro interpro.tsv

smips interpro.tsv

# smips2cassis.pl 40 interpro.tsv.anchor_genes.csv <annotation> <genome>


cat interpro.tsv.anchor_genes.csv | grep -v "^#" | cut -f1 > anchor_genes.txt

cat $Genes | grep 'mRNA' | sed 's/ID=//g' | sed "s/;.*//g" | awk '{ print $9 "\t" $1 "\t" $4 "\t" $5 "\t" $7}' > cassis.tsv

for GeneID in $(cat anchor_genes.txt | tail -n+2); do
  echo $GeneID
  mkdir -p $GeneID
  cassis \
    --annotation cassis.tsv \
    --genome $Assembly \
    --anchor $GeneID \
    --dir out/$GeneID \
    --mismatches 0 \
    -v \
    --prediction \
    --num-cpus 40 \
    | tee 2>&1 out/$GeneID/${GeneID}_log.txt
done
```

Results for the Fusarin cluster were summarised:

```bash
for File in $(ls out/g2085.t1/*/CLUSTER/*/*/fimo/fimo.txt); do
  Gene=$(echo $File | cut -f2 -d '/')
  Region=$(echo $File | cut -f5 -d '/')
  UpDown=$(echo $File | cut -f6 -d '/')
  Motif=$(cat $File | tail -n+2 |cut -f1 | sort | uniq)
  printf "${Gene}\t${Region}\t${UpDown}\t${Motif}\n"
done | sort -k4
```

This did not include the Tri5 cluster, so these genes were run seperately:

```bash
for GeneID in g3123.t1 g3124.t1 g3125.t1 g3126.t1 g3127.t1 g3128.t1 g3129.t1 g3130.t1 g3131.t1 g3132.t1 g3133.t1 g3134.t1 g3135.t1; do
  echo $GeneID
  mkdir -p tri_cluster_out/$GeneID
  cassis \
    --annotation cassis.tsv \
    --genome $Assembly \
    --anchor $GeneID \
    --dir tri_cluster_out/$GeneID \
    --mismatches 0 \
    -v \
    --prediction \
    --num-cpus 40 \
    | tee 2>&1 tri_cluster_out/$GeneID/${GeneID}_log.txt
done
```

Results were summarised:

```bash
for File in $(ls tri_cluster_out/*/*/CLUSTER/*/*/fimo/fimo.txt); do
  Gene=$(echo $File | cut -f2 -d '/')
  Region=$(echo $File | cut -f5 -d '/')
  UpDown=$(echo $File | cut -f6 -d '/')
  Motif=$(cat $File | tail -n+2 |cut -f1 | sort | uniq)
  printf "${Gene}\t${Region}\t${UpDown}\t${Motif}\n"
done | sort -k4
```

```
g3124.t1        g3124.t1_to_g3132.t1    +6_-6   ABGCTGTACTC
g3125.t1        g3124.t1_to_g3132.t1    +6_-6   ABGCTGTACTC
g3127.t1        g3124.t1_to_g3132.t1    +4_-8   ABGCTGTACTC
g3128.t1        g3124.t1_to_g3132.t1    +4_-8   ABGCTGTACTC
g3129.t1        g3124.t1_to_g3132.t1    +3_-9   ABGCTGTACTC
g3132.t1        g3124.t1_to_g3132.t1    +0_-12  ABGCTGTACTC
g3123.t1        g3116.t1_to_g3131.t1    +2_-4   CACTTTAVCCC
g3126.t1        g3118.t1_to_g3134.t1    +6_-1   GRNAGGCCTDR
g3130.t1        g3118.t1_to_g3134.t1    +3_-4   GRNAGGCCTDR
g3130.t1        g3124.t1_to_g3134.t1    downstream_border_+3_-4 GRNAGGCCTDR
g3131.t1        g3118.t1_to_g3134.t1    +2_-5   GRNAGGCCTDR
g3131.t1        g3118.t1_to_g3135.t1    upstream_border_+2_-5   GRNAGGCCTDR
g3131.t1        g3124.t1_to_g3134.t1    downstream_border_+2_-5 GRNAGGCCTDR
g3133.t1        g3118.t1_to_g3135.t1    upstream_border_+0_-7   GRNAGGCCTDR
g3134.t1        g3118.t1_to_g3135.t1    upstream_border_+0_-7   GRNAGGCCTDR
g3130.t1        g3124.t1_to_g3134.t1    upstream_border_+4_-1   TGGTHGGCCTA
g3131.t1        g3118.t1_to_g3135.t1    downstream_border_+3_-2 TGGTHGGCCTA
g3131.t1        g3124.t1_to_g3134.t1    upstream_border_+3_-2   TGGTHGGCCTA
g3131.t1        g3124.t1_to_g3135.t1    +3_-2   TGGTHGGCCTA
g3133.t1        g3118.t1_to_g3135.t1    downstream_border_+1_-4 TGGTHGGCCTA
g3133.t1        g3124.t1_to_g3135.t1    +1_-4   TGGTHGGCCTA
g3134.t1        g3118.t1_to_g3135.t1    downstream_border_+1_-4 TGGTHGGCCTA
g3134.t1        g3124.t1_to_g3135.t1    +1_-4   TGGTHGGCCTA
g3135.t1        g3124.t1_to_g3135.t1    +0_-5   TGGTHGGCCTA
```

## Confirmation of results

MEME was run on the Tri5 region to replicate the results:

```bash
Promoters=$(ls out/g*.t1/PROMOTERS/all_promoter_sequences.fasta | head -n1)
OutDir=tri_meme_out
mkdir -p $OutDir/meme
TriPromoters=tri_meme_out/tri_promoters.fa
cat $Promoters | sed -n '/^>g3124.t1/,/^>g3136.t1/p' | grep -v "^$" |head -n -1 > $TriPromoters

meme $TriPromoters -dna -mod anr -nmotifs 1 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc $OutDir/meme

cat tri_meme_out/meme/meme.txt | grep -C2 'regular expression'

#---
# Running FIMO
#---
Thresh=0.00006 #CASSIS (trained for secmet clusters)
# Thresh=1e-4 #Default
fimo -thresh $Thresh -oc tri_meme_out/fimo tri_meme_out/meme/meme.html $Promoters
echo "Promoters containing motif:"
cat tri_meme_out/fimo/fimo.txt | cut -f3 | sort | uniq | wc -l
echo "Total number found in these promoters:"
cat tri_meme_out/fimo/fimo.txt | cut -f3 | sort | wc -l
echo "this covers the following genes:"
cat tri_meme_out/fimo/fimo.txt | cut -f3 | sort | uniq | grep -o -P "g.*?\.t\d" > tri_meme_out/fimo/gene_containing_motifs.txt
cat tri_meme_out/fimo/gene_containing_motifs.txt | wc -l

AnnotTab=$(ls /oldhpc/home/groups/harrisonlab/project_files/fusarium_venenatum/gene_pred/annotation/F.venenatum/WT/WT_annotation_ncbi_expression.tsv)
for GeneID in $(cat tri_meme_out/fimo/gene_containing_motifs.txt); do
cat $AnnotTab | grep  "^${GeneID}"
done > tri_meme_out/fimo/gene_containing_motifs_annots.tsv

cat tri_meme_out/fimo/gene_containing_motifs_annots.tsv | grep -i 'DEG:' | wc -l
cat $AnnotTab | grep -i 'DEG:' | wc -l
# 183 of 2199 DEGs contain this motif.


#---
# Running MAST
#---

mast tri_meme_out/meme/meme.xml $TriPromoters -oc $OutDir/mast -comp
ls tri_meme_out/mast/mast.txt
```


Attempt to identify gapped motifs in Tri cluster:

```bash
OutDir=tri_glam_out
mkdir -p $OutDir/glam
TriPromoters=tri_meme_out/tri_promoters.fa
glam2 n -O $OutDir/glam -2 -n 10000 $TriPromoters
ls tri_glam_out/glam/glam2.txt
```


Attempt to identify gapped motifs Tri6 & Tri10 genes:

```bash
WorkDir=~/tmp/cassis_Fv
mkdir $WorkDir
cd $WorkDir
Promoters=$(ls out/g*.t1/PROMOTERS/all_promoter_sequences.fasta | head -n1)
OutDir=tri_glam_out
mkdir -p $OutDir/glam

GeneList="g3128.t1 g3130.t1"
)
SubPromoters=$OutDir/Tri6_Tri10_promoters.fa

for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>$GeneID/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

glam2 n -O $OutDir/glam_Tri6_Tri10 -2 -n 10000 $SubPromoters
ls $OutDir/glam_Tri6_Tri10/glam2.txt
```


### Fusarin Cluster

Attempt to identify common elements in the Fusarin gene cluster.


```bash
screen -a
ssh compute02

WorkDir=~/tmp/cassis_Fv
mkdir $WorkDir
cd $WorkDir
Promoters=$(ls out/g*.t1/PROMOTERS/all_promoter_sequences.fasta | head -n1)
OutDir=fusarin_meme_out
mkdir -p $OutDir/meme

GeneList="g2081.t1 g2082.t1 g2083.t1 g2079.t1 g2077.t1 g2084.t1 g2076.t1 g2085.t1 g2078.t1 g2080.t1 g2060.t1"
# GeneList=$(for num in $(seq 2032 2095); do printf "g${num}.t1 "; done)

ListLen=$(echo $GeneList | grep -o 'g' | wc -l)
SubPromoters=fusarin_meme_out/fusarin_promoters_${ListLen}.fa

for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>$GeneID/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

meme $SubPromoters -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc $OutDir/meme_${ListLen}

cat $OutDir/meme_${ListLen}/meme.txt | grep -C2 'regular expression' | grep -v "^--$" | grep -v "^$"

# From g3130.t1 g1528.t1"
# Motif GGMCATRCAAAG MEME-1 regular expression
# GG[AC]CAT[AG]CAAAG

#---
# Running FIMO
#---
Thresh=0.00006 #CASSIS (trained for secmet clusters)
# Thresh=1e-4 #Default
fimo -thresh $Thresh -oc $OutDir/fimo_${ListLen} $OutDir/meme_${ListLen}/meme.html fusarin_meme_out/fusarin_promoters_64.fa
for Motif in $(cat $OutDir/fimo_${ListLen}/fimo.tsv | grep 'contig' | cut -f1 | sort | uniq); do
  echo $Motif
  echo "Promoters containing motif:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | wc -l
  echo "Total number found in these promoters:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | wc -l
  echo "this covers the following genes:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | grep -o -P "g.*?\.t\d" > $OutDir/fimo_${ListLen}/gene_containing_motifs.txt
  cat $OutDir/fimo_${ListLen}/gene_containing_motifs.txt | wc -l
done

# AnnotTab=$(ls /oldhpc/home/groups/harrisonlab/project_files/fusarium_venenatum/gene_pred/annotation/F.venenatum/WT/WT_annotation_ncbi_expression.tsv)
# for GeneID in $(cat tri_meme_out/fimo/gene_containing_motifs.txt); do
# cat $AnnotTab | grep  "^${GeneID}"
# done > tri_meme_out/fimo/gene_containing_motifs_annots.tsv
#
# cat tri_meme_out/fimo/gene_containing_motifs_annots.tsv | grep -i 'DEG:' | wc -l
# cat $AnnotTab | grep -i 'DEG:' | wc -l
#


#---
# Running MAST
#---

mast $OutDir/meme_${ListLen}/meme.xml $SubPromoters -oc $OutDir/mast
mast $OutDir/meme_${ListLen}/meme.xml fusarin_meme_out/fusarin_promoters_64.fa -oc $OutDir/mast_all_cluster
ls $OutDir/mast/mast.txt
```

```
--------------------------------------------------------------------------------
        Motif ATCCTTCR MEME-1 regular expression
--------------------------------------------------------------------------------
ATC[CT]T[TA]C[AG]

ATCCTTCR
Promoters containing motif:
7
Total number found in these promoters:
12
this covers the following genes:
11
```


```
--------------------------------------------------------------------------------
	Motif YYGAAGGAT MEME-1 regular expression
--------------------------------------------------------------------------------
[TC][CT]GAAGGAT
--------------------------------------------------------------------------------
	Motif CGCMRCACCMG MEME-2 regular expression
--------------------------------------------------------------------------------
CG[CG][AC][AG][CGT]ACC[AC]G
--------------------------------------------------------------------------------
	Motif GCGCCGCCRCC MEME-3 regular expression
--------------------------------------------------------------------------------
GCGCCGCC[AG]CC
--------------------------------------------------------------------------------
	Motif KGCAADTGTCAG MEME-4 regular expression
--------------------------------------------------------------------------------
[GT][GC][CG]A[AG][TAG]TG[TG]CAG
--------------------------------------------------------------------------------
	Motif CGGADGGTTTYA MEME-5 regular expression
--------------------------------------------------------------------------------
CGGA[GAT][GT]GTTT[CT]A

CGCMRCACCMG
Promoters containing motif:
5
Total number found in these promoters:
9
this covers the following genes:
8
CGGADGGTTTYA
Promoters containing motif:
3
Total number found in these promoters:
7
this covers the following genes:
4
GCGCCGCCRCC
Promoters containing motif:
2
Total number found in these promoters:
2
this covers the following genes:
4
KGCAADTGTCAG
Promoters containing motif:
6
Total number found in these promoters:
7
this covers the following genes:
9
YYGAAGGAT
Promoters containing motif:
6
Total number found in these promoters:
15
this covers the following genes:
9
```
<!--


### Transcription factors

Attempt to identify common elements in similarly regulated transcription factors.


```bash
screen -a
ssh compute02


WorkDir=~/tmp/cassis_Fv
mkdir $WorkDir
cd $WorkDir
Promoters=$(ls out/g*.t1/PROMOTERS/all_promoter_sequences.fasta | head -n1)
OutDir=tri_meme_TF_out
mkdir -p $OutDir/meme

# GeneList="g10.t1 g1967.t1 g1528.t1 g3130.t1 g10338.t1"
# GeneList="g1704.t1 g12379.t1 g12427.t1 g1490.t1 g2437.t1 g6597.t1 g7360.t1 g1888.t1 g9268.t1 g8487.t1 g9623.t1 g2182.t1 g7676.t1 g9017.t1 g6382.t1 g10234.t1 g6306.t1 g1541.t1 g820.t1 g10390.t1 g8825.t1 g499.t1 g7497.t1 g4653.t1 g5510.t1 g8341.t1 g1290.t1 g9924.t1 g9275.t1 g10558.t1 g1407.t1 g10.t1 g1967.t1 g1528.t1 g3130.t1 g10338.t1 g4714.t1 g6017.t1 g320.t1 g12195.t1 g3653.t1 g3508.t1 g2035.t1 g10056.t1 g1037.t1 g11177.t1 g5257.t1 g10694.t1 g2100.t1 g5166.t1"
GeneList="g3130.t1 g1528.t1"
ListLen=$(echo $GeneList | grep -o 'g' | wc -l)
SubPromoters=tri_meme_TF_out/TF_promoters_${ListLen}.fa

for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>$GeneID/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

meme $SubPromoters -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc $OutDir/meme_${ListLen}

cat $OutDir/meme_${ListLen}/meme.txt | grep -C2 'regular expression'

# From g3130.t1 g1528.t1"
# Motif GGMCATRCAAAG MEME-1 regular expression
# GG[AC]CAT[AG]CAAAG

#---
# Running FIMO
#---
Thresh=0.00006 #CASSIS (trained for secmet clusters)
# Thresh=1e-4 #Default
fimo -thresh $Thresh -oc $OutDir/fimo_${ListLen} $OutDir/meme_${ListLen}/meme.html $OutDir/TF_promoters_50.fa
echo "Promoters containing motif:"
cat $OutDir/fimo_${ListLen}/fimo.tsv | cut -f3 | sort | uniq | wc -l
echo "Total number found in these promoters:"
cat $OutDir/fimo_${ListLen}/fimo.tsv | cut -f3 | sort | wc -l
echo "this covers the following genes:"
cat $OutDir/fimo_${ListLen}/fimo.tsv | cut -f3 | sort | uniq | grep -o -P "g.*?\.t\d" > $OutDir/fimo_${ListLen}/gene_containing_motifs.txt
cat $OutDir/fimo_${ListLen}/gene_containing_motifs.txt | wc -l

# AnnotTab=$(ls /oldhpc/home/groups/harrisonlab/project_files/fusarium_venenatum/gene_pred/annotation/F.venenatum/WT/WT_annotation_ncbi_expression.tsv)
# for GeneID in $(cat tri_meme_out/fimo/gene_containing_motifs.txt); do
# cat $AnnotTab | grep  "^${GeneID}"
# done > tri_meme_out/fimo/gene_containing_motifs_annots.tsv
#
# cat tri_meme_out/fimo/gene_containing_motifs_annots.tsv | grep -i 'DEG:' | wc -l
# cat $AnnotTab | grep -i 'DEG:' | wc -l
#


#---
# Running MAST
#---

mast $OutDir/meme_${ListLen}/meme.xml $SubPromoters -oc $OutDir/mast
ls $OutDir/mast/mast.txt
```

FIMO identified 21 occurences in 15 genes from 50 highly expressed TFs.

The analysis was repeated using these 15 genes.

```bash
screen -a
ssh compute02


WorkDir=~/tmp/cassis_Fv
mkdir $WorkDir
cd $WorkDir
Promoters=$(ls out/g*.t1/PROMOTERS/all_promoter_sequences.fasta | head -n1)
OutDir=tri_meme_TF_out
mkdir -p $OutDir/meme

ListLen=$(cat tri_meme_TF_out/fimo_2/fimo.tsv | grep 'contig' | cut -f3 | sort | uniq | wc -l)
SubPromoters=tri_meme_TF_out/TF_promoters_${ListLen}.fa

for GeneID in $(cat tri_meme_TF_out/fimo_2/fimo.tsv | grep 'contig' | cut -f3 | sort | uniq); do
cat $Promoters | sed -n "/^>$GeneID/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

meme $SubPromoters -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc $OutDir/meme_${ListLen}

cat $OutDir/meme_${ListLen}/meme.txt | grep -C2 'regular expression'

# From g3130.t1 g1528.t1"
# Motif GGMCATRCAAAG MEME-1 regular expression
# GG[AC]CAT[AG]CAAAG

#---
# Running MAST
#---

mast $OutDir/meme_${ListLen}/meme.xml $SubPromoters -oc $OutDir/mast
ls $OutDir/mast/mast.txt

#---
# Running FIMO
#---
Thresh=0.00006 #CASSIS (trained for secmet clusters)
# Thresh=1e-4 #Default
fimo -thresh $Thresh -oc $OutDir/fimo_${ListLen} $OutDir/meme_${ListLen}/meme.html $OutDir/TF_promoters_50.fa
for Motif in $(cat $OutDir/fimo_${ListLen}/fimo.tsv | grep 'contig' | cut -f1 | sort | uniq); do
  echo $Motif
  echo "Promoters containing motif:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | wc -l
  echo "Total number found in these promoters:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | wc -l
  echo "this covers the following genes:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | grep -o -P "g.*?\.t\d" > $OutDir/fimo_${ListLen}/gene_containing_motifs.txt
  cat $OutDir/fimo_${ListLen}/gene_containing_motifs.txt | wc -l
done
```

In set of 50 genes
```
GGCRGTCKNGGC
Promoters containing motif:
14
Total number found in these promoters:
18
this covers the following genes:
15
TCTTTCTTCTTY
Promoters containing motif:
22
Total number found in these promoters:
56
this covers the following genes:
24
```


Attempt to identify gapped motifs in TF similarly expressed genes:

```bash
WorkDir=~/tmp/cassis_Fv
mkdir $WorkDir
cd $WorkDir
Promoters=$(ls out/g*.t1/PROMOTERS/all_promoter_sequences.fasta | head -n1)
OutDir=tri_glam_TF_out
mkdir -p $OutDir/glam

# GeneList="g1704.t1 g12379.t1 g12427.t1 g1490.t1 g2437.t1 g6597.t1 g7360.t1 g1888.t1 g9268.t1 g8487.t1 g9623.t1 g2182.t1 g7676.t1 g9017.t1 g6382.t1 g10234.t1 g6306.t1 g1541.t1 g820.t1 g10390.t1 g8825.t1 g499.t1 g7497.t1 g4653.t1 g5510.t1 g8341.t1 g1290.t1 g9924.t1 g9275.t1 g10558.t1 g1407.t1 g10.t1 g1967.t1 g1528.t1 g3130.t1 g10338.t1 g4714.t1 g6017.t1 g320.t1 g12195.t1 g3653.t1 g3508.t1 g2035.t1 g10056.t1 g1037.t1 g11177.t1 g5257.t1 g10694.t1 g2100.t1 g5166.t1"
# GeneList="g10.t1 g1967.t1 g1528.t1 g3130.t1 g10338.t1"
GeneList="g3130.t1 g1528.t1"
ListLen=$(echo $GeneList | grep -o 'g' | wc -l)
SubPromoters=$OutDir/TF_promoters_${ListLen}.fa

for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>$GeneID/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

glam2 n -O $OutDir/glam_${ListLen} -2 -n 10000 $SubPromoters
ls $OutDir/glam_${ListLen}/glam2.txt
```
 -->




# Analysis using UTR gene set:


### Tri genes

Attempt to identify common elements in similarly regulated transcription factors.


```bash
screen -a
ssh compute02


WorkDir=~/tmp/cassis_Fv
mkdir $WorkDir
cd $WorkDir
ProjDir=$(ls -d /oldhpc/home/groups/harrisonlab/project_files/fusarium_venenatum)
Promoters=$(ls $ProjDir/analysis/promoters/F.venenatum/WT_UTR/WT_promoters.fa)
OutDir=tri_meme_UTR_out
mkdir -p $OutDir/meme

GeneList="g3123 g3124 g3125 g3126 g3127 g3128 g3129 g3130 g3131 g3132 g3133 g3134 g3135"
ListLen=$(echo $GeneList | grep -o 'g' | wc -l)
SubPromoters=$OutDir/TF_promoters_${ListLen}.fa

Old2NewTable=$(ls $ProjDir/gene_pred/braker/*/*_UTR/*/gene_conversion/gene_intersects.gff)
for GeneID in $GeneList; do
  NewGeneID=$(cat $Old2NewTable | grep "ID=${GeneID};.contig" | cut -f18 | sed 's/ID=//g' | tr -d ';')
cat $Promoters | sed -n "/^>${NewGeneID}_/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

meme $SubPromoters -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc $OutDir/meme_${ListLen}

cat $OutDir/meme_${ListLen}/meme.txt | grep -C2 'regular expression'


#---
# Running MAST
#---

mast $OutDir/meme_${ListLen}/meme.xml $SubPromoters -oc $OutDir/mast
ls $OutDir/mast/mast.txt

#---
# Running FIMO
#---
Thresh=0.00006 #CASSIS (trained for secmet clusters)
# Thresh=1e-4 #Default
fimo -thresh $Thresh -oc $OutDir/fimo_${ListLen} $OutDir/meme_${ListLen}/meme.html $OutDir/TF_promoters_${ListLen}.fa
for Motif in $(cat $OutDir/fimo_${ListLen}/fimo.tsv | grep 'TSS' | cut -f1 | sort | uniq); do
  echo $Motif
  echo "Promoters containing motif:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | wc -l
  echo "Total number found in these promoters:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | wc -l
  echo "this covers the following genes:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | grep -o -P "g.*?\.t\d" > $OutDir/fimo_${ListLen}/gene_containing_motifs.txt
  cat $OutDir/fimo_${ListLen}/gene_containing_motifs.txt | wc -l
done

ls -d $PWD/$OutDir
```

```
CTTTCMMCYGCM
Promoters containing motif:
8
Total number found in these promoters:
29
this covers the following genes:
0
GGKRGGCCTDR
Promoters containing motif:
10
Total number found in these promoters:
26
this covers the following genes:
0
RATRAWCCATCT
Promoters containing motif:
7
Total number found in these promoters:
16
this covers the following genes:
0
TAAWCCAKGGCY
Promoters containing motif:
4
Total number found in these promoters:
10
this covers the following genes:
0
TGRGGCTGTWMT
Promoters containing motif:
6
Total number found in these promoters:
8
this covers the following genes:
0
```


### Transcription factors

Attempt to identify common elements in similarly regulated transcription factors.


```bash
screen -a
ssh compute02


WorkDir=~/tmp/cassis_Fv
mkdir $WorkDir
cd $WorkDir
# Promoters=$(ls out/g*.t1/PROMOTERS/all_promoter_sequences.fasta | head -n1)
ProjDir=$(ls -d /oldhpc/home/groups/harrisonlab/project_files/fusarium_venenatum)
Promoters=$(ls $ProjDir/analysis/promoters/F.venenatum/WT_UTR/WT_promoters.fa)
OutDir=tri_meme_UTR_TF_out
mkdir -p $OutDir/meme

# GeneList="g10.t1 g1967.t1 g1528.t1 g3130.t1 g10338.t1"

GeneList="g1704 g12379 g12427 g1490 g2437 g6597 g7360 g1888 g9268 g8487 g9623 g2182 g7676 g9017 g6382 g10234 g6306 g1541 g820 g10390 g8825 g499 g7497 g4653 g5510 g8341 g1290 g9924 g9275 g10558 g1407 g10 g1967 g1528 g3130 g10338 g4714 g6017 g320 g12195 g3653 g3508 g2035 g10056 g1037 g11177 g5257 g10694 g2100 g5166"
ListLen=$(echo $GeneList | grep -o 'g' | wc -l)
SubPromoters=$OutDir/TF_promoters_${ListLen}.fa

Old2NewTable=$(ls $ProjDir/gene_pred/braker/*/*_UTR/*/gene_conversion/gene_intersects.gff)
for GeneID in $GeneList; do
  NewGeneID=$(cat $Old2NewTable | grep "ID=${GeneID};.contig" | cut -f18 | sed 's/ID=//g' | tr -d ';')
  # printf "${GeneID} $NewGeneID\n"
cat $Promoters | sed -n "/^>${NewGeneID}_/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

GeneList="g3130 g1528"
ListLen=$(echo $GeneList | grep -o 'g' | wc -l)
SubPromoters=$OutDir/TF_promoters_${ListLen}.fa

Old2NewTable=$(ls $ProjDir/gene_pred/braker/*/*_UTR/*/gene_conversion/gene_intersects.gff)
for GeneID in $GeneList; do
  NewGeneID=$(cat $Old2NewTable | grep "ID=${GeneID};.contig" | cut -f18 | sed 's/ID=//g' | tr -d ';')
  # echo "${GeneID} - ${NewGeneID}"
cat $Promoters | sed -n "/^>${NewGeneID}_/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

meme $SubPromoters -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc $OutDir/meme_${ListLen}

cat $OutDir/meme_${ListLen}/meme.txt | grep -C2 'regular expression'

# From g3130.t1 g1528.t1"
# Motif GGMCATRCAAAG MEME-1 regular expression
# GG[AC]CAT[AG]CAAAG

#---
# Running MAST
#---

mast $OutDir/meme_${ListLen}/meme.xml $SubPromoters -oc $OutDir/mast
ls $OutDir/mast/mast.txt

#---
# Running FIMO
#---
# Thresh=0.00006 #CASSIS (trained for secmet clusters)
Thresh=1e-4 #Default
fimo -thresh $Thresh -oc $OutDir/fimo_${ListLen} $OutDir/meme_${ListLen}/meme.html $OutDir/TF_promoters_50.fa
for Motif in $(cat $OutDir/fimo_${ListLen}/fimo.tsv | grep 'TSS' | cut -f1 | sort | uniq); do
  echo $Motif
  echo "Promoters containing motif:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | wc -l
  echo "Total number found in these promoters:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | wc -l
  echo "this covers the following genes:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | grep -o -P "g\d+" > $OutDir/fimo_${ListLen}/${Motif}_gene_containing_motifs_newnames.txt
  cat $OutDir/fimo_${ListLen}/${Motif}_gene_containing_motifs_newnames.txt | wc -l
  Old2NewTable=$(ls $ProjDir/gene_pred/braker/*/*_UTR/*/gene_conversion/gene_intersects.gff)
  for NewName in $(cat $OutDir/fimo_${ListLen}/${Motif}_gene_containing_motifs_newnames.txt); do
    OldGeneID=$(cat $Old2NewTable | grep "ID=${NewName};$" | cut -f9 | sed 's/ID=//g' | tr -d ';')
    printf "${OldGeneID}\t$NewName\n"
  done > $OutDir/fimo_${ListLen}/${Motif}_gene_containing_motifs.txt
done
```
```
CAMATGCA
Promoters containing motif:
9
Total number found in these promoters:
15
this covers the following genes:
9
GGMCATRCAAAG
Promoters containing motif:
10
Total number found in these promoters:
16
this covers the following genes:
10
KWCAGAAANCA
Promoters containing motif:
11
Total number found in these promoters:
15
this covers the following genes:
11
TTTTCTCSAWS
Promoters containing motif:
12
Total number found in these promoters:
18
this covers the following genes:
12

```


```
[A/T][A/T]CAAAG binding motif of HMG (high mobility group) TFs:
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0015199
```


Attempt to identify common elements in similarly regulated transcription factors.


```bash
screen -a
ssh compute02


WorkDir=~/tmp/cassis_Fv
mkdir $WorkDir
cd $WorkDir
# Promoters=$(ls out/g*.t1/PROMOTERS/all_promoter_sequences.fasta | head -n1)
ProjDir=$(ls -d /oldhpc/home/groups/harrisonlab/project_files/fusarium_venenatum)
Promoters=$(ls $ProjDir/analysis/promoters/F.venenatum/WT_UTR/WT_promoters.fa)
OutDir=meme_UTR_TF_out
mkdir -p $OutDir/meme

# GeneList="g10.t1 g1967.t1 g1528.t1 g3130.t1 g10338.t1"

GeneList="g1704 g12379 g12427 g1490 g2437 g6597 g7360 g1888 g9268 g8487 g9623 g2182 g7676 g9017 g6382 g10234 g6306 g1541 g820 g10390 g8825 g499 g7497 g4653 g5510 g8341 g1290 g9924 g9275 g10558 g1407 g10 g1967 g1528 g3130 g10338 g4714 g6017 g320 g12195 g3653 g3508 g2035 g10056 g1037 g11177 g5257 g10694 g2100 g5166"
ListLen=$(echo $GeneList | grep -o 'g' | wc -l)
SubPromoters=$OutDir/TF_promoters_${ListLen}.fa

Old2NewTable=$(ls $ProjDir/gene_pred/braker/*/*_UTR/*/gene_conversion/gene_intersects.gff)
for GeneID in $GeneList; do
  NewGeneID=$(cat $Old2NewTable | grep "ID=${GeneID};.contig" | cut -f18 | sed 's/ID=//g' | tr -d ';')
  # printf "${GeneID} $NewGeneID\n"
cat $Promoters | sed -n "/^>${NewGeneID}_/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

GeneList="g6017 g3653 g1407 g1528 g2035 g3130 g8825 g12379 g9275 g4714"
ListLen=$(echo $GeneList | grep -o 'g' | wc -l)
SubPromoters=$OutDir/TF_promoters_${ListLen}.fa

Old2NewTable=$(ls $ProjDir/gene_pred/braker/*/*_UTR/*/gene_conversion/gene_intersects.gff)
for GeneID in $GeneList; do
  NewGeneID=$(cat $Old2NewTable | grep "ID=${GeneID};.contig" | cut -f18 | sed 's/ID=//g' | tr -d ';')
  # echo "${GeneID} - ${NewGeneID}"
cat $Promoters | sed -n "/^>${NewGeneID}_/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

meme $SubPromoters -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc $OutDir/meme_${ListLen}

cat $OutDir/meme_${ListLen}/meme.txt | grep -C2 'regular expression'

#---
# Running MAST
#---

mast $OutDir/meme_${ListLen}/meme.xml $SubPromoters -oc $OutDir/mast_${ListLen}
ls $OutDir/mast/mast.txt

mast $OutDir/meme_${ListLen}/meme.xml $OutDir/TF_promoters_50.fa -oc $OutDir/mast_${ListLen}_vs_50

#---
# Running FIMO
#---
# Thresh=0.00006 #CASSIS (trained for secmet clusters)
Thresh=1e-4 #Default
fimo -thresh $Thresh -oc $OutDir/fimo_${ListLen} $OutDir/meme_${ListLen}/meme.html $OutDir/TF_promoters_50.fa
for Motif in $(cat $OutDir/fimo_${ListLen}/fimo.tsv | grep 'TSS' | cut -f1 | sort | uniq); do
  echo $Motif
  echo "Promoters containing motif:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | wc -l
  echo "Total number found in these promoters:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | wc -l
  echo "this covers the following genes:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | grep -o -P "g\d+" > $OutDir/fimo_${ListLen}/${Motif}_gene_containing_motifs_newnames.txt
  cat $OutDir/fimo_${ListLen}/${Motif}_gene_containing_motifs_newnames.txt | wc -l
  Old2NewTable=$(ls $ProjDir/gene_pred/braker/*/*_UTR/*/gene_conversion/gene_intersects.gff)
  for NewName in $(cat $OutDir/fimo_${ListLen}/${Motif}_gene_containing_motifs_newnames.txt); do
    OldGeneID=$(cat $Old2NewTable | grep "ID=${NewName};$" | cut -f9 | sed 's/ID=//g' | tr -d ';')
    printf "${OldGeneID}\t$NewName\n"
  done > $OutDir/fimo_${ListLen}/${Motif}_gene_containing_motifs.txt
done
```



### Fusarin genes

Attempt to identify common elements in similarly regulated transcription factors.


```bash
screen -a
ssh compute02


WorkDir=~/tmp/cassis_Fv
mkdir $WorkDir
cd $WorkDir
ProjDir=$(ls -d /oldhpc/home/groups/harrisonlab/project_files/fusarium_venenatum)
Promoters=$(ls $ProjDir/analysis/promoters/F.venenatum/WT_UTR/WT_promoters.fa)
OutDir=fusarin_meme_UTR_out
mkdir -p $OutDir/meme

GeneList=$(for num in $(seq 2032 2095); do printf "g${num} "; done)
ListLen=$(echo $GeneList | grep -o 'g' | wc -l)
SubPromoters=$OutDir/TF_promoters_${ListLen}.fa

Old2NewTable=$(ls $ProjDir/gene_pred/braker/*/*_UTR/*/gene_conversion/gene_intersects.gff)
for GeneID in $GeneList; do
  NewGeneID=$(cat $Old2NewTable | grep "ID=${GeneID};.contig" | cut -f18 | sed 's/ID=//g' | tr -d ';')
  # printf "${GeneID} $NewGeneID\n"
cat $Promoters | sed -n "/^>${NewGeneID}_/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

GeneList="g2081 g2082 g2083 g2079 g2077 g2084 g2076 g2085 g2078 g2080 g2060"
ListLen=$(echo $GeneList | grep -o 'g' | wc -l)
SubPromoters=$OutDir/TF_promoters_${ListLen}.fa

Old2NewTable=$(ls $ProjDir/gene_pred/braker/*/*_UTR/*/gene_conversion/gene_intersects.gff)
for GeneID in $GeneList; do
  NewGeneID=$(cat $Old2NewTable | grep "ID=${GeneID};.contig" | cut -f18 | sed 's/ID=//g' | tr -d ';')
  # echo $GeneID
cat $Promoters | sed -n "/^>${NewGeneID}_/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

meme $SubPromoters -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc $OutDir/meme_${ListLen}

cat $OutDir/meme_${ListLen}/meme.txt | grep -C2 'regular expression'


#---
# Running MAST
#---

mast $OutDir/meme_${ListLen}/meme.xml $SubPromoters -oc $OutDir/mast
ls $OutDir/mast/mast.txt

#---
# Running FIMO
#---
Thresh=0.00006 #CASSIS (trained for secmet clusters)
# Thresh=1e-4 #Default
fimo -thresh $Thresh -oc $OutDir/fimo_${ListLen} $OutDir/meme_${ListLen}/meme.html $OutDir/TF_promoters_64.fa
for Motif in $(cat $OutDir/fimo_${ListLen}/fimo.tsv | grep 'TSS' | cut -f1 | sort | uniq); do
  echo $Motif
  echo "Promoters containing motif:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | wc -l
  echo "Total number found in these promoters:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | wc -l
  echo "this covers the following genes:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | grep -o -P "g.*?\.t\d" > $OutDir/fimo_${ListLen}/gene_containing_motifs.txt
  cat $OutDir/fimo_${ListLen}/gene_containing_motifs.txt | wc -l
done
```

```
ACTGTCARACTC
Promoters containing motif:
9
Total number found in these promoters:
15
this covers the following genes:
0
ARATCCTTCGR
Promoters containing motif:
17
Total number found in these promoters:
37
this covers the following genes:
0
ARWAGAARRAKA
Promoters containing motif:
17
Total number found in these promoters:
25
this covers the following genes:
0
GATCTTATCTYA
Promoters containing motif:
15
Total number found in these promoters:
21
this covers the following genes:
0
TGCABWCRCCC
Promoters containing motif:
10
Total number found in these promoters:
12
this covers the following genes:
0
```
