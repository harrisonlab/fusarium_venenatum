# Generating an TSV file with sequencing information

```bash
#Antismash output correction
cat analysis/secondary_metabolites/antismash/F.venenatum/WT_minion_VP/WT_antismash_results_secmet_genes.tsv | sed 's/;//p' | sed 's/;.*//p' | sed 's/Kin.*//p' > analysis/secondary_metabolites/antismash/F.venenatum/WT_minion_VP/WT_antismash_results_secmet_genes_corrected.tsv


for GeneGff in $(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.gff3); do
Strain=WT_minion
Organism=F.venenatum
Assembly=$(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
TFs=$(ls analysis/transcription_factors/F.venenatum/WT_minion/WT_minion_TF_domains.tsv)
InterPro=$(ls gene_pred/interproscan/F.venenatum/WT_minion/WT_minion_interproscan.tsv)
Antismash=$(ls analysis/secondary_metabolites/antismash/F.venenatum/WT_minion_VP/WT_antismash_results_secmet_genes_corrected.tsv)
#Smurf=$(ls analysis/secondary_metabolites/smurf/F.venenatum/WT_minion/WT_minion_smurf_secmet_genes.tsv) # I added cassis genes manually
SwissProt=$(ls gene_pred/swissprot/F.venenatum/WT_minion/swissprot_vJun2020_tophit_parsed.tbl)
OutDir=analysis/annotation_tables_VP/$Organism/$Strain
mkdir -p $OutDir
GeneFasta=$(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta)
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Annotation_tables
$ProgDir/build_annot_v1.py --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --TFs $TFs --InterPro $InterPro --Antismash $Antismash --Swissprot $SwissProt > $OutDir/"$Strain"_noDEGs_gene_table.tsv
done
```


## TRI5 promoters

MEME was run on the Tri5 region to replicate the results:

```bash
# Identify sequence motifs in the TRI5 cluster promoters
meme tri5_cluster/tri_promoters.fa -dna -mod anr -nmotifs 3 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc tri5_cluster/meme
cat tri5_cluster/meme/meme.txt | grep -C2 'regular expression'

# Scan motifs in all promoter sequences
Thresh=0.00006 #CASSIS (trained for secmet clusters), thresh=1e-4 #Default
fimo -thresh $Thresh -oc tri5_cluster/fimo tri5_cluster/meme/meme.html all_promoter_sequences.fasta
for Motif in $(cat tri5_cluster/fimo/fimo.tsv | grep 'contig' | cut -f1 | sort | uniq); do
echo $Motif
echo "Promoters containing motif:"
cat tri5_cluster/fimo/fimo.tsv | cut -f3 | sort | uniq | wc -l
echo "Total number found in these promoters:"
cat tri5_cluster/fimo/fimo.tsv | cut -f3 | sort | wc -l
echo "this covers the following genes:"
cat tri5_cluster/fimo/fimo.tsv | cut -f3 | sort | uniq | grep -o -P "g.*?\.t\d" > tri5_cluster/fimo/gene_containing_motifs.txt
cat tri5_cluster/fimo/gene_containing_motifs.txt | wc -l
```
```
GRDAGGCCTRA
Promoters containing motif:
1068
Total number found in these promoters:
1172
this covers the following genes:
1368
```
```bash
mast tri5_cluster/meme/meme.html tri5_cluster/tri_promoters.fa -oc tri5_cluster/mast -comp

# Ame needs an older version of meme and python 2.7
conda create --name meme_py27 python=2.7
conda install 'meme=5.0.2' 'icu=58.2'
# Control tri5 promoters
ame --control tri5_cluster/tri_promoters.fa --oc tri5_cluster/ame all_promoter_sequences.fasta tri5_cluster/meme/meme.html
# No control
ame --oc tri5_cluster/ame_nocontrol all_promoter_sequences.fasta tri5_cluster/meme/meme.html

# Identify gapped motifs. This needs same ame installation
glam2 n -O tri5_cluster/glam -2 -n 10000 tri5_cluster/tri_promoters.fa
less tri5_cluster/glam/glam2.txt





# Identify gapped motifs Tri6 & Tri10 genes:
Promoters=$(ls all_promoter_sequences.fasta)
GeneList="g6429.t1 g6432.t1"
SubPromoters=Tri6_Tri10_promoters.fa
for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>$GeneID/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters
glam2 n -O tri5_cluster/glam_Tri6_Tri10 -2 -n 10000 $SubPromoters
ls $OutDir/glam_Tri6_Tri10/glam2.txt


# AnnotTab=$(ls /projects/oldhome/groups/harrisonlab/project_files/fusarium_venenatum/gene_pred/annotation/F.venenatum/WT/WT_annotation_ncbi_expression.tsv)
# for GeneID in $(cat tri_meme_out/fimo/gene_containing_motifs.txt); do
# cat $AnnotTab | grep  "^${GeneID}"
# done > tri_meme_out/fimo/gene_containing_motifs_annots.tsv

# cat tri_meme_out/fimo/gene_containing_motifs_annots.tsv | grep -i 'DEG:' | wc -l
# cat $AnnotTab | grep -i 'DEG:' | wc -l
# 183 of 2199 DEGs contain this motif.
```

Analysis of a subset of genes

```bash
# Subset 1
Promoters=$(ls all_promoter_sequences.fasta)
OutDir=tri5_subsets
mkdir -p $OutDir/TF_g4107/meme

GeneList="g6426.t1 g6428.t1 g6429.t1 g6430.t1 g6431.t1 g6433.t1 g6435.t1"
# GeneList=$(for num in $(seq 2032 2095); do printf "g${num}.t1 "; done)

for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>$GeneID/,/^>/p" | grep -v "^$" | head -n -1
done > $OutDir/TF_g4107/promoters.fasta

meme $OutDir/TF_g4107/promoters.fasta -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc $OutDir/TF_g4107/meme

# Subset 2
Promoters=$(ls all_promoter_sequences.fasta)
OutDir=tri5_subsets
mkdir -p $OutDir/TF_g4106/meme

GeneList="g6426.t1 g6427.t1 g6428.t1 g6429.t1 g6430.t1 g6431.t1 g6434.t1 g6435.t1 g6436.t1 g6437.t1 g8466.t1"
# GeneList=$(for num in $(seq 2032 2095); do printf "g${num}.t1 "; done)

for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>$GeneID/,/^>/p" | grep -v "^$" | head -n -1
done > $OutDir/TF_g4106/promoters.fasta

meme $OutDir/TF_g4106/promoters.fasta -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc $OutDir/TF_g4106/meme

# Subset 3
Promoters=$(ls all_promoter_sequences.fasta)
OutDir=tri5_subsets
mkdir -p $OutDir/TF_g6432/meme

GeneList="g6426.t1 g6430.t1 g6431.t1 g6433.t g6434.t1 g6435.t1 g6436.t1"
# GeneList=$(for num in $(seq 2032 2095); do printf "g${num}.t1 "; done)

for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>$GeneID/,/^>/p" | grep -v "^$" | head -n -1
done > $OutDir/TF_g6432/promoters.fasta

meme $OutDir/TF_g6432/promoters.fasta -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc $OutDir/TF_g6432/meme
# This gave a different motif with the lowest P-val

# Subset 4. Same as before but larger max width
Promoters=$(ls all_promoter_sequences.fasta)
OutDir=tri5_subsets
mkdir -p $OutDir/TF_g6432_maxw12/meme

GeneList="g6426.t1 g6430.t1 g6431.t1 g6433.t g6434.t1 g6435.t1 g6436.t1"
# GeneList=$(for num in $(seq 2032 2095); do printf "g${num}.t1 "; done)

for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>$GeneID/,/^>/p" | grep -v "^$" | head -n -1
done > $OutDir/TF_g6432_maxw12/promoters.fasta

meme $OutDir/TF_g6432_maxw12/promoters.fasta -dna -mod anr -nmotifs 5 -minw 6 -maxw 20 -revcomp -evt 1.0e+005 -oc $OutDir/TF_g6432_maxw12/meme

```


## Fusarin Cluster

Attempt to identify common elements in the Fusarin gene cluster.


```bash
# Meme
Promoters=$(ls all_promoter_sequences.fasta)
OutDir=fusarin_cluster
mkdir -p $OutDir/meme

GeneList="g12320.t1 g12337.t1 g12338.t1 g12339.t1 g12340.t1 g12341.t1 g12342.t1 g12343.t1 g12344.t1 g12345.t1 g12346.t1"
# GeneList=$(for num in $(seq 2032 2095); do printf "g${num}.t1 "; done)

ListLen=$(echo $GeneList | grep -o 'g' | wc -l)
SubPromoters=fusarin_cluster/fusarin_promoters_${ListLen}.fa

for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>$GeneID/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

meme $SubPromoters -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc $OutDir/meme_${ListLen}
cat $OutDir/meme_${ListLen}/meme.txt | grep -C2 'regular expression' | grep -v "^--$" | grep -v "^$"

#Not significant

# Fimo
Thresh=0.00006 #CASSIS (trained for secmet clusters)
# Thresh=1e-4 #Default
fimo -thresh $Thresh -oc $OutDir/fimo_${ListLen} $OutDir/meme_${ListLen}/meme.html fusarin_cluster/fusarin_promoters_11.fa
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
```
ATCCTTCR
Promoters containing motif:
5
Total number found in these promoters:
10
this covers the following genes:
8
```

```bash
fimo -thresh $Thresh -oc $OutDir/fimo_all $OutDir/meme_${ListLen}/meme.html all_promoter_sequences.fasta
for Motif in $(cat $OutDir/fimo_all/fimo.tsv | grep 'contig' | cut -f1 | sort | uniq); do
  echo $Motif
  echo "Promoters containing motif:"
  cat $OutDir/fimo_all/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | wc -l
  echo "Total number found in these promoters:"
  cat $OutDir/fimo_all/fimo.tsv | grep "^$Motif" | cut -f3 | sort | wc -l
  echo "this covers the following genes:"
  cat $OutDir/fimo_all/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | grep -o -P "g.*?\.t\d" > $OutDir/fimo_all/gene_containing_motifs.txt
  cat $OutDir/fimo_all/gene_containing_motifs.txt | wc -l
done
```
```
ATCCTTCR
Promoters containing motif:
1105
Total number found in these promoters:
1226
this covers the following genes:
1456
```

```bash
# Mast
mast $OutDir/meme_${ListLen}/meme.xml $SubPromoters -oc $OutDir/mast --comp
# This can be done with other secmet genes but no black space between genes is needed.
```







# AnnotTab=$(ls /projects/oldhome/groups/harrisonlab/project_files/fusarium_venenatum/gene_pred/annotation/F.venenatum/WT/WT_annotation_ncbi_expression.tsv)
# for GeneID in $(cat tri_meme_out/fimo/gene_containing_motifs.txt); do
# cat $AnnotTab | grep  "^${GeneID}"
# done > tri_meme_out/fimo/gene_containing_motifs_annots.tsv
#
# cat tri_meme_out/fimo/gene_containing_motifs_annots.tsv | grep -i 'DEG:' | wc -l
# cat $AnnotTab | grep -i 'DEG:' | wc -l


```
AA promoters . More than me, probably gene promoters are separated and shorter 
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



### Transcription factors

Attempt to identify common elements in similarly regulated transcription factors.


```bash
Promoters=$(ls all_promoter_sequences.fasta)
OutDir=tri_meme_TF_out
mkdir -p $OutDir/meme

GeneList="g10479 g4106 g4107.t1 g4134 g6128 g6432 g66.t1 g5413 g6132"
ListLen=$(echo $GeneList | grep -o 'g' | wc -l)
SubPromoters=tri_meme_TF_out/TF_promoters_${ListLen}.fa

for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>$GeneID/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

meme $SubPromoters -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc $OutDir/meme_${ListLen}

cat $OutDir/meme_${ListLen}/meme.txt | grep -C2 'regular expression' | grep -v "^--$" | grep -v "^$"
```
```
--------------------------------------------------------------------------------
        Motif CMTCAAS MEME-1 regular expression
--------------------------------------------------------------------------------
C[AC]TCAA[GC]
--------------------------------------------------------------------------------
        Motif ACCAADMCAGG MEME-2 regular expression
--------------------------------------------------------------------------------
A[CT]C[AC]A[AGT][CA]C[AG]GG
```
Not significant
```bash
# Fimo
Thresh=0.00006 #CASSIS (trained for secmet clusters)
# Thresh=1e-4 #Default
fimo -thresh $Thresh -oc $OutDir/fimo_${ListLen} $OutDir/meme_${ListLen}/meme.html $OutDir/TF_promoters_9.fa
for Motif in $(cat tri_meme_TF_out/fimo_9/fimo.tsv | grep 'contig' | cut -f1 | sort | uniq); do
echo $Motif
echo "Promoters containing motif:"
cat tri_meme_TF_out/fimo_9/fimo.tsv | cut -f3 | sort | uniq | wc -l
echo "Total number found in these promoters:"
cat tri_meme_TF_out/fimo_9/fimo.tsv | cut -f3 | sort | wc -l
echo "this covers the following genes:"
cat tri_meme_TF_out/fimo_9/fimo.tsv | cut -f3 | sort | uniq | grep -o -P "g.*?\.t\d" > tri_meme_TF_out/fimo_9/gene_containing_motifs.txt
cat tri_meme_TF_out/fimo_9/gene_containing_motifs.txt | wc -l
done
```
```
ACCAADMCAGG
Promoters containing motif:
10
Total number found in these promoters:
17
this covers the following genes:
8
```
