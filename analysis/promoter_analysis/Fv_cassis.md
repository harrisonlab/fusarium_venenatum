

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


Results for each cluster were summarised:

```bash
cd ~/tmp/cassis_Fv
AnnotTab=$(ls /projects/oldhome/groups/harrisonlab/project_files/fusarium_venenatum/gene_pred/annotation/F.venenatum/WT/WT_annotation_ncbi_expression.tsv)
for Results in $(ls out/*/*_log.txt | grep -v 'g6640'); do
Anchor=$(echo $Results | cut -f2 -d '/')
SortID=$(echo $Anchor | tr -d /g/ | cut -f1 -d '.')
Cluster=$(cat $AnnotTab | cut -f1,7 | grep "$Anchor" | cut -f2)
if [[ $Cluster == "" ]]; then
  Cluster="NA"
fi
if grep 'No cluster prediction' $Results; then
printf "${SortID}\t${Anchor}\t${Cluster}\tNA\tNA\n"
else
Best=$(cat $Results | grep -A2 '(7) Computing CLUSTER PREDICTIONS' | tail -n1 | sed -r "s&^\s+&&g" | cut -f1 -d ' ')
Fimo=$(ls out/$Anchor/$Anchor/fimo/$Best/fimo.txt)
Motif=$(cat $Fimo | head -n2 | tail -n1 | cut -f1)
printf "${SortID}\t${Anchor}\t${Cluster}\t${Best}\t${Motif}\n"
fi
done | grep -v ':(' | sort -n -k1 | cut -f2-
```

```
g717.t1	NA	+2_-3	BTCAGCTCRY
g755.t1	SecMet_cluster_2	+2_-3	CGCYAYGGCCCG
g1114.t1	SecMet_cluster_6	+10_-4	ACGATGTCCA
g1181.t1	SecMet_cluster_7	+3_-7	GYCCRCCMVYMC
g1968.t1	NA	NA	NA
g1976.t1	NA	+2_-4	TCAAMGAACC
g2085.t1	SecMet_cluster_10	+1_-6	YYGAAGGATC
g2143.t1	SecMet_cluster_11	+15_-9	CCCASTAGCYTM
g2729.t1	NA	+0_-3	AGACTTTRMA
g2817.t1	SecMet_cluster_13	+1_-2	CACAARCTGCC
g2902.t1	SecMet_cluster_14	+0_-10	CGCCTGGMATTG
g3649.t1	SecMet_cluster_17	+4_-6	YTCGGANKCCG
g3842.t1	SecMet_cluster_18	+1_-6	TGGGGTGTCWTG
g3971.t1	SecMet_cluster_19	+0_-3	GARCSGCATMTT
g4399.t1	SecMet_cluster_20	+2_-2	SGGACCGAAGG
g4406.t1	SecMet_cluster_20	+4_-0	ACCCKRCGGAC
g4464.t1	SecMet_cluster_21	+3_-0	TCAAGGTCSTC
g4783.t1	NA	+3_-0	YGSRMGSGAG
g4921.t1	SecMet_cluster_23	NA	NA
g5005.t1	NA	+1_-8	CATCSGCYATGG
g5605.t1	SecMet_cluster_24	+0_-15	AAAACTACTATT
g6008.t1	NA	+1_-2	MCAGATGCATC
g6374.t1	NA	+0_-5	AKVGYRGGGAT
g6376.t1	NA	+15_-7	AAAATTGCAGTG
g6417.t1	SecMet_cluster_26	+0_-3	CYKTGTCTTAC
g6459.t1	NA	+14_-1	CAACCCCTSCT
g6648.t1	SecMet_cluster_29	+0_-4	CGYGAAGGAWAA
g7068.t1	SecMet_cluster_31	+5_-0	CBTKATTCAA
g7393.t1	NA	+0_-8	CGACAAWKCSAA
g7496.t1	NA	+2_-3	CAAWGCCKAC
g7509.t1	SecMet_cluster_32	+2_-2	TGWTATTGACA
g7513.t1	SecMet_cluster_32	+0_-4	TGWTATTGACA
g8047.t1	SecMet_cluster_34	+2_-11	ACTGGGACTRGG
g8510.t1	NA	+0_-4	TYGGRGCTYCGG
g8877.t1	SecMet_cluster_37	+5_-0	GCCGCGGA
g9021.t1	NA	+3_-1	RGTGKGGKWGKC
g9129.t1	NA	+5_-6	TGGATCTTCRAG
g9373.t1	SecMet_cluster_40	+0_-6	TCAAGKCATCSK
g9515.t1	SecMet_cluster_41	+1_-2	KCYGTTGRTC
g9634.t1	NA	+7_-0	SCATGGC
g10233.t1	NA	+9_-9	CGYSRCTGGCG
g10406.t1	NA	+2_-9	AAGCCWCTCGMT
g10578.t1	NA	+2_-7	TACBTACCTAC
g11163.t1	SecMet_cluster_47	NA	NA
g11164.t1	NA	+6_-11	ACACTMCCCTM
g11190.t1	NA	+3_-2	CCAGAWGCCATG
g12093.t1	SecMet_cluster_49	+14_-12	TAARTAGYATAA
g12106.t1	SecMet_cluster_49	+6_-14	TTATGCTATTTA
g12148.t1	SecMet_cluster_50	+3_-3	YGGRGAAWTGRG
g12225.t1	NA	+2_-1	CMCCAAGTCTA
g12540.t1	SecMet_cluster_51	+0_-4	GGMCCCGBCCCC
g12555.t1	NA	+1_-2	TGABCGTCGCC
```

Cassis results were copied to the project directory

```bash
  WorkDir=/newhome/armita/tmp/cassis_Fv
  OutDir=/home/groups/harrisonlab/project_files/fusarium_venenatum/analysis/secondary_metabolites/cassis
  mkdir $OutDir

  cp -r $WorkDir/out $OutDir/.
  cp -r $WorkDir/*anchor_gene* $OutDir/.
  cp -r $WorkDir/cassis.tsv $OutDir/.
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

### Cassis on all genes in genome


```bash
# screen -a
# ssh compute02
# srun -p medium --mem=4gb --pty bash

conda activate meme-v4
OldProjDir=$(ls -d /projects/oldhome/groups/harrisonlab/project_files/fusarium_venenatum)

# WorkDir=~/tmp/cassis_Fv/all_genes
WorkDir=$OldProjDir/analysis/promoters/cassis/all_genes
mkdir -p $WorkDir
cd $WorkDir

AnnotTab=$(ls $OldProjDir/gene_pred/annotation/F.venenatum/WT/WT_annotation_ncbi_expression.tsv)
Assembly=$(ls $OldProjDir/repeat_masked/F.venenatum/WT/illumina_assembly_ncbi/WT_contigs_softmasked_repeatmasker_TPSI_appended.fa)
Genes=$(ls $OldProjDir/gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.gff3)
Interpro=$(ls $OldProjDir/gene_pred/interproscan/F.venenatum/WT/WT_interproscan.tsv)

cat $Genes | grep 'mRNA' | sed 's/ID=//g' | sed "s/;.*//g" | awk '{ print $9 "\t" $1 "\t" $4 "\t" $5 "\t" $7}' > cassis.tsv

for Cluster in $(cat $AnnotTab | cut -f7 | grep 'SecMet_cluster' | sort -n -k3 -t'_' | uniq); do
echo $Cluster
mkdir $WorkDir/$Cluster
cat $AnnotTab | cut -f1,7 | grep -w "$Cluster" | cut -f1 | grep '.t1' > $WorkDir/$Cluster/headers.txt
for GeneID in $(cat $WorkDir/$Cluster/headers.txt); do
  echo $GeneID
  mkdir -p $WorkDir/$Cluster/$GeneID
  ProgDir=/projects/oldhome/armita/git_repos/emr_repos/scripts/fusarium_venenatum/analysis/promoter_analysis
  CassisTSV=cassis.tsv
  OutDir=$WorkDir/$Cluster
  Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
  while [ $Jobs -gt 60 ]; do
    sleep 5m
    printf "."
    Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
  done		
  printf "\n"
  sbatch $ProgDir/cassis_SLURM.sh $Assembly $CassisTSV $GeneID $OutDir
  # cassis \
  #   --annotation cassis.tsv \
  #   --genome $Assembly \
  #   --anchor $GeneID \
  #   --dir $WorkDir/$Cluster/$GeneID \
  #   --mismatches 0 \
  #   -v \
  #   --prediction \
  #   --num-cpus 30 \
  #   | tee 2>&1 $WorkDir/$Cluster/$GeneID/${GeneID}_log.txt
done
done
```

```bash
ProjDir=$(ls -d /projects/oldhome/groups/harrisonlab/project_files/fusarium_venenatum)
cd $ProjDir
for Cluster in $(ls -d analysis/promoters/cassis/all_genes/SecMet_cluster_* | rev | cut -f1 -d '/' | rev | sort -n -k3 -t'_'); do
ClusterDir=$(ls -d analysis/promoters/cassis/all_genes/${Cluster})
echo ""
for Results in $(ls $ClusterDir/*/*_log.txt); do
Anchor=$(echo $Results | rev | cut -f2 -d '/' | rev)
# if [[ $Cluster == "" ]]; then
# Cluster="NA"
# fi
if $(grep -q 'No cluster prediction' $Results); then
printf "${Cluster}\t${Anchor}\tNA\tNA\n"
elif grep 'Computing CLUSTER PREDICTIONS' $Results; then
Best=$(cat $Results | grep -A2 '(7) Computing CLUSTER PREDICTIONS' | tail -n1 | sed -r "s&^\s+&&g" | cut -f1 -d ' ')
Fimo=$(ls $ClusterDir/$Anchor/$Anchor/fimo/$Best/fimo.txt)
Motif=$(cat $Fimo | head -n2 | tail -n1 | cut -f1)
printf "${Cluster}\t${Anchor}\t${Best}\t${Motif}\n"
else
  printf "${Cluster}\t${Anchor}\tNA\tNA\n"
fi
done | grep -v 'CLUSTER PREDICTIONS' | grep -v ':('
done > analysis/promoters/cassis/all_genes/cassis_summary.tsv
```

SecMet genes with transcription factor annotations were identified

```bash
ProjDir=$(ls -d /projects/oldhome/groups/harrisonlab/project_files/fusarium_venenatum)
AnnotTab=$(ls $ProjDir/gene_pred/annotation/F.venenatum/WT/WT_annotation_ncbi_expression.tsv)
cat $AnnotTab | grep 'SecMet_cluster' | cut -f1,7,14 | grep -v  "\s$"
```

Differentially expressed SecMet genes were identified:

```bash
ProjDir=$(ls -d /projects/oldhome/groups/harrisonlab/project_files/fusarium_venenatum)
AnnotTab=$(ls $ProjDir/gene_pred/annotation/F.venenatum/WT/WT_annotation_ncbi_expression.tsv)
cat $AnnotTab | grep 'SecMet_cluster' | cut -f1,7,49 | grep -v  "\s$"
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

Attempt to identify common elements in similarly regulated Tri cluster genes.


```bash
screen -a
# ssh compute02
srun -p medium --mem=4gb --pty bash


WorkDir=~/tmp/cassis_Fv
mkdir $WorkDir
cd $WorkDir
ProjDir=$(ls -d /projects/oldhome/groups/harrisonlab/project_files/fusarium_venenatum)
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
# ssh compute02
# srun -p medium --mem=4gb --pty bash
qlogin -pe smp=4


# Promoters=$(ls out/g*.t1/PROMOTERS/all_promoter_sequences.fasta | head -n1)
# ProjDir=$(ls -d /oldhpc/home/groups/harrisonlab/project_files/fusarium_venenatum)
ProjDir=$(ls -d /home/groups/harrisonlab/project_files/fusarium_venenatum)
Promoters=$(ls $ProjDir/analysis/promoters/F.venenatum/WT_UTR/WT_promoters.fa)

# WorkDir=~/tmp/cassis_Fv
WorkDir=$ProjDir
mkdir $WorkDir
cd $WorkDir

# OutDir=tri_meme_UTR_TF_out
OutDir=analysis/promtoers/TF/meme_UTR_TF_out
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

meme $SubPromoters -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc $OutDir/meme_${ListLen}

cat $OutDir/meme_${ListLen}/meme.txt | grep -C2 'regular expression'

# Poor motif - just T's

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

mast $OutDir/meme_${ListLen}/meme.xml $SubPromoters -oc $OutDir/meme_${ListLen}/mast_${ListLen}
mast $OutDir/meme_${ListLen}/meme.xml $OutDir/TF_promoters_50.fa -oc $OutDir/meme_${ListLen}/mast_50

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
