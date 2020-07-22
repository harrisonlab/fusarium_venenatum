# Generating an TSV file with sequencing information

```bash
#Antismash output correction
cat /projects/fusarium_venenatum/analysis/secondary_metabolites/antismash/F.venenatum/WT_minion_vAG/WT_antismash_results_secmet_genes.tsv | sed 's/;//p' | sed 's/;.*//p' | sed 's/Kin.*//p' > WT_antismash_results_secmet_genes_corrected.tsv


for GeneGff in $(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.gff3); do
    Strain=WT_minion
    Organism=F.venenatum
    Assembly=$(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
    TFs=$(ls analysis/transcription_factors/F.venenatum/WT_minion/WT_minion_TF_domains.tsv)
    InterPro=$(ls gene_pred/interproscan/F.venenatum/WT_minion/WT_minion_interproscan.tsv)
    Antismash=$(ls analysis/secondary_metabolites/antismash/F.venenatum/WT_minion_vAG/WT_antismash_results_secmet_genes_corrected.tsv)
    SwissProt=$(ls gene_pred/swissprot/F.venenatum/WT_minion/swissprot_vJun2020_tophit_parsed.tbl)
    OutDir=analysis/annotation_tables/$Organism/$Strain
    mkdir -p $OutDir
    GeneFasta=$(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta)
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Annotation_tables
    $ProgDir/build_annot_v1.py --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --TFs $TFs --InterPro $InterPro --Antismash $Antismash --Swissprot $SwissProt > $OutDir/"$Strain"_noDEGs_gene_table.tsv
done
```

# Cassis on all genes 

```bash

# meme version 4 required
conda activate meme_v4

ProjDir=$(ls -d /projects/fusarium_venenatum)

WorkDir=$ProjDir/analysis/promoters/cassis/all_genes
mkdir -p $WorkDir
cd $WorkDir

AnnotTab=$(ls $ProjDir/analysis/annotation_tables/F.venenatum/WT_minion/WT_minion_noDEGs_gene_table.tsv)
Assembly=$(ls $ProjDir/repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/medaka_contigs_softmasked_repeatmasker_TPSI_appended.fa)
Genes=$(ls $ProjDir/gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.gff3)
Interpro=$(ls $ProjDir/gene_pred/interproscan/F.venenatum/WT_minion/WT_minion_interproscan.tsv)

cat $Genes | grep 'mRNA' | sed 's/ID=//g' | sed "s/;.*//g" | awk '{ print $9 "\t" $1 "\t" $4 "\t" $5 "\t" $7}' > cassis.tsv

CassisTSV=cassis.tsv

for Cluster in $(cat $AnnotTab | cut -f13 | grep 'contig' | sort -n -k3 -t'_' | sed 's/;.*//p' | uniq); do # Extract cluster names
echo $Cluster
mkdir $WorkDir/$Cluster
cat $AnnotTab | cut -f1,13 | grep -w "$Cluster" | cut -f1 | grep '.t1' > $WorkDir/$Cluster/headers.txt # Create a headers file per cluster
for GeneID in $(cat $WorkDir/$Cluster/headers.txt); do
echo $GeneID
mkdir -p $WorkDir/$Cluster/$GeneID
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Promoter_analysis
OutDir=$WorkDir/$Cluster
Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
while [ $Jobs -gt 60 ]; do
sleep 5m
printf "."
Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
done
printf "\n"
sbatch $ProgDir/cassis.sh $Assembly $CassisTSV $GeneID $OutDir
done
done
```

```bash
ProjDir=$(ls -d /projects/fusarium_venenatum)
cd $ProjDir
for Cluster in $(ls -d analysis/promoters/cassis/all_genes3/contig* | rev | cut -f1 -d '/' | rev | sort -n -k3 -t'_'); do
ClusterDir=$(ls -d analysis/promoters/cassis/all_genes3/${Cluster})
echo ""
for Results in $(ls $ClusterDir/*/*_log.txt); do
Anchor=$(echo $Results | rev | cut -f2 -d '/' | rev)
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
done > analysis/promoters/cassis/all_genes3/cassis_summary.tsv
```
```
Tri5 cluster
contig_2_Cluster_10     g6427.t1        +8_-0   GRDAGGCCTRA
contig_2_Cluster_10     g6428.t1        +7_-1   GRDAGGCCTRA
contig_2_Cluster_10     g6429.t1        +6_-2   GRDAGGCCTRA
contig_2_Cluster_10     g6430.t1        +6_-2   GRDAGGCCTRA
contig_2_Cluster_10     g6431.t1        +5_-3   GRDAGGCCTRA
contig_2_Cluster_10     g6432.t1        +4_-1   GRKAGGCCTDR
contig_2_Cluster_10     g6433.t1        +3_-0   YYAGGCCTCYC
contig_2_Cluster_10     g6434.t1        +2_-1   YYAGGCCTCYC
contig_2_Cluster_10     g6435.t1        +1_-2   YYAGGCCTCYC
contig_2_Cluster_10     g6436.t1        +0_-3   YYAGGCCTCYC
```

cat analysis/annotation_tables/F.venenatum/WT_minion/WT_minion_noDEGs_gene_table.tsv | grep 'Cluster' | cut -f1,11,13 | grep -v  "\s$"



```bash
# Extract best motifs meme format

ProjDir=$(ls -d /projects/fusarium_venenatum)
cd $ProjDir
for Cluster in $(ls -d analysis/promoters/cassis/all_genes3/contig* | rev | cut -f1 -d '/' | rev | sort -n -k3 -t'_'); do
ClusterDir=$(ls -d analysis/promoters/cassis/all_genes3/${Cluster})
for Results in $(ls $ClusterDir/*/*_log.txt); do
Anchor=$(echo $Results | rev | cut -f2 -d '/' | rev)
if $(grep -q 'No cluster prediction' $Results); then
continue
elif grep 'Computing CLUSTER PREDICTIONS' $Results; then
Best=$(cat $Results | grep -A2 '(7) Computing CLUSTER PREDICTIONS' | tail -n1 | sed -r "s&^\s+&&g" | cut -f1 -d ' ')
Best_ed=$(echo $Best | tr -d '+' | tr -d '-')
# Fimo=$(ls $ClusterDir/$Anchor/$Anchor/fimo/$Best/fimo.txt)
# Motif=$(cat $Fimo | head -n2 | tail -n1 | cut -f1)
MotifFile=$(ls $ClusterDir/$Anchor/$Anchor/meme/$Best/meme.txt)
echo '--------------------------------------------------------------------------------'
cat $MotifFile | awk '/position-specific probability matrix/,/expression/' | head -n-4 | sed "s/MEME-1/${Anchor}_${Best_ed}_MEME-1/g"
echo ""
else
continue
fi
done | grep -v 'CLUSTER PREDICTIONS' | grep -v ':('
done > analysis/promoters/cassis/all_genes/best_motifs_meme.txt



ProjDir=$(ls -d /projects/fusarium_venenatum)
cd $ProjDir
printf "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n" > analysis/promoters/cassis/all_genes/best_motifs_meme_for_mast.txt
for Cluster in $(ls -d analysis/promoters/cassis/all_genes/contig* | rev | cut -f1 -d '/' | rev | sort -n -k3 -t'_'); do
ClusterDir=$(ls -d analysis/promoters/cassis/all_genes/${Cluster})
for Results in $(ls $ClusterDir/*/*_log.txt); do
Anchor=$(echo $Results | rev | cut -f2 -d '/' | rev)
if $(grep -q 'No cluster prediction' $Results); then
continue
elif grep 'Computing CLUSTER PREDICTIONS' $Results; then
Best=$(cat $Results | grep -A2 '(7) Computing CLUSTER PREDICTIONS' | tail -n1 | sed -r "s&^\s+&&g" | cut -f1 -d ' ')
Best_ed=$(echo $Best | tr -d '+' | tr -d '-')
MotifFile=$(ls $ClusterDir/$Anchor/$Anchor/meme/$Best/meme.txt)
#echo '--------------------------------------------------------------------------------'
cat $MotifFile | awk '/position-specific probability matrix/,/expression/' | head -n-4 | sed "s/MEME-1/${Anchor}_${Best_ed}_MEME-1/g"
echo ""
else
continue
fi
done
done | grep -v 'CLUSTER PREDICTIONS' | grep -v ':(' >> analysis/promoters/cassis/all_genes/best_motifs_meme_for_mast.txt

```

```
Please, enter number of files to read (must be > 0):
1
Please, enter number of top significant motifs (must be > 0 and <= 50):
10
Please, enter number of best matches (must be > 0 and <= 50):
10
Please, select a cutoff for similarity (>= 0.5, >= 0.6, >= 0.7, >= 0.75, >= 0.8, >= 0.85, >= 0.9):
0.8
Please, enter number of threads (must be >= 1):
4
Maximum number of threads available on your machine is 24.
This is the maximum number of threads can be allocated to run this program.

Please, enter input file's location (full path, for example, C:\MyDocuments\ for Windows and /home/MyFolder/ for Linux):
/projects/fusarium_venenatum/analysis/promoters/cassis/all_genes3/

Enter input file names and formats (for example: 1). See the user manual for each format:

(1) TRANSFAC
(2) TRANSFAC-like
(3) PSSM
(4) Jaspar
(5) MEME output
(6) Consensus sequence
(7) Sequence Alignment
(8) Matrices (Horizonal)
(9) Matrices (Vertical)
(10) Unspecified

Please, enter file name (in text format .txt, name without spaces):
best_motifs_meme.txt
Please, enter file format:
5

Enter database name by selecting a number in the list below:

(1) Jaspar 2016 (All)
(2) Jaspar 2016 (Fungi)
(3) Jaspar 2016 (Insects)
(4) Jaspar 2016 (Nematodes)
(5) Jaspar 2016 (Plants)
(6) Jaspar 2016 (Urochordates)
(7) Jaspar 2016 (Vertebrates)
(8) Transfac (Free version)
(9) UniPROBE (Human)
(10) UniPROBE (Mouse)
(11) UniPROBE (Parasite)
(12) UniPROBE (Worm)
(13) UniPROBE (Yeast)
(14) None

2
Would you like to generate motif tree? (Y or N):
Y
Would you like to combine similar motifs? (Y or N):
Y
Please, enter an output file type (Global-Only, All):
All
Please, enter an output file format (Text, HTML, PDF, All):
All
please, enter output file's location (full path, for example, C:\MyDocuments\ for Windows and /home/MyFolder/ for Linux):
/projects/fusarium_venenatum/analysis/promoters/cassis/
```

Motif similarity 

```bash
cd /projects/fusarium_venenatum
mkdir analysis/promoters/cassis/all_genes3/motif_similarity
printf "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n" > analysis/promoters/cassis/all_genes3/motif_similarity/all_motifs.txt
for Cluster in $(ls -d analysis/promoters/cassis/all_genes3/contig* | rev | cut -f1 -d '/' | rev | sort -n -k3 -t'_'); do
ClusterDir=$(ls -d analysis/promoters/cassis/all_genes3/${Cluster})
for Results in $(ls $ClusterDir/*/*_log.txt); do
Anchor=$(echo $Results | rev | cut -f2 -d '/' | rev)
if $(grep -q 'No cluster prediction' $Results); then
continue
elif grep 'Computing CLUSTER PREDICTIONS' $Results; then
Best=$(cat $Results | grep -A2 '(7) Computing CLUSTER PREDICTIONS' | tail -n1 | sed -r "s&^\s+&&g" | cut -f1 -d ' ')
Best_ed=$(echo $Best | tr -d '+' | tr -d '-')
MotifFile=$(ls $ClusterDir/$Anchor/$Anchor/meme/$Best/meme.txt)
cat $MotifFile | grep "^MOTIF" | grep -v "MOTIF DIAGRAM" | sed "s/MOTIF/MOTIF ${Cluster}_${Anchor}_${Best}/g" | sed "s/MEME-1\\s//g"
cat $MotifFile | awk '/letter-probability matrix/,/---/' | head -n-1
echo ""
else
continue
fi
done
done | grep -v 'CLUSTER PREDICTIONS' | grep -v ':(' >> analysis/promoters/cassis/all_genes3/motif_similarity/all_motifs.txt
```

```bash
AllMotifs=$(ls analysis/promoters/cassis/all_genes/motif_similarity/all_motifs.txt)
# Create a file with the following lines
printf "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n" > analysis/promoters/cassis/all_genes/out/best_motifs_meme.txt
for Cluster in $(ls -d analysis/promoters/cassis/all_genes/contig* | rev | cut -f1 -d '/' | rev | sort -n -k3 -t'_'); do
ClusterDir=$(ls -d analysis/promoters/cassis/all_genes/${Cluster})
echo "$Cluster"
OutDir=analysis/promoters/cassis/all_genes/motif_similarity/$Cluster
mkdir $OutDir
cat $AllMotifs | awk "/${Cluster}/,/^$/{next}{print}" > $OutDir/all_motifs_excl_${Cluster}.txt
for Results in $(ls $ClusterDir/*/*_log.txt); do
Anchor=$(echo $Results | rev | cut -f2 -d '/' | rev)
if $(grep -q 'No cluster prediction' $Results); then
continue
elif grep 'Computing CLUSTER PREDICTIONS' $Results; then
Best=$(cat $Results | grep -A2 '(7) Computing CLUSTER PREDICTIONS' | tail -n1 | sed -r "s&^\s+&&g" | cut -f1 -d ' ')
MotifFile=$(ls $ClusterDir/$Anchor/$Anchor/meme/$Best/meme.txt)
mkdir $OutDir/$Anchor
tomtom -oc $OutDir/$Anchor/tomtom $MotifFile $OutDir/all_motifs_excl_${Cluster}.txt
echo ""
else
continue
fi
done
done

for Cluster in $(ls -d analysis/promoters/cassis/all_genes/motif_similarity/contig_* | rev | cut -f1 -d '/' | rev | sort -n -k3 -t'_'); do
ClusterDir=$(ls -d analysis/promoters/cassis/all_genes/motif_similarity/${Cluster})
for TomTom in $(ls $ClusterDir/*/tomtom/tomtom.txt); do
Anchor=$(echo $TomTom | rev | cut -f3 -d '/' | rev)
cat $TomTom | cut -f1,2,4,9,10 | grep 'Cluster' | grep -v '#' | sed "s/_g.*_-[0-9]*//g" | sort | uniq | sed "s/^/$Cluster\t$Anchor\t/g"
done
done > analysis/promoters/cassis/all_genes/motif_similarity/tomtom_hits.tsv
cat analysis/promoters/cassis/all_genes/motif_similarity/tomtom_hits.tsv | grep 'e-' > analysis/promoters/cassis/all_genes/motif_similarity/tomtom_hits_high_score.tsv

for Cluster in $(ls -d analysis/promoters/cassis/all_genes/motif_similarity/contig_* | rev | cut -f1 -d '/' | rev | sort -n -k3 -t'_'); do
ClusterDir=$(ls -d analysis/promoters/cassis/all_genes/motif_similarity/${Cluster})
for TomTom in $(ls $ClusterDir/*/tomtom/tomtom.txt); do
Anchor=$(echo $TomTom | rev | cut -f3 -d '/' | rev)
cat $TomTom | grep 'Cluster' | grep -v '#' | sed "s/_g.*_-[0-9]*//g" | sort | uniq | sed "s/^/$Cluster\t$Anchor\t/g"
done
done > analysis/promoters/cassis/all_genes/motif_similarity/tomtom_hits_full.tsv
```



MAST
```bash
cd analysis/promoters/mast
/data/scratch/gomeza/prog/smips WT_minion_interproscan.tsv 
cat WT_minion_interproscan.tsv.anchor_genes.csv |  grep -v "^#" | cut -f1 > anchor_genes.txt

Promoters=$(ls analysis/promoters/cassis/all_genes/contig_2_Cluster_8/g6184.t1/PROMOTERS/all_promoter_sequences.fasta)
for AnchorGene in $(cat analysis/promoters/mast/anchor_genes.txt); do
cat $Promoters | awk "/$AnchorGene/,/^$/" | head -n-1
done > analysis/promoters/anchor_genes.fa
```




mast analysis/promoters/cassis/all_genes/motif_similarity/all_motifs.txt analysis/promoters/mast/anchor_genes.fa -oc analysis/promoters/mast -comp

analysis/promoters/cassis/all_genes/best_motifs_meme_for_mast.txt 

-----------------

Extract promotor regions:
Promotor regions upstream of all genes were extracted:

```bash
for GeneGff in $(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.gff3); do
Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
Assembly=$(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
OutDir=analysis/meme/promotor_regions/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Promoter_analysis
$ProgDir/extract_promotor.pl --fasta $Assembly --gff $GeneGff --prefix $OutDir/${Strain} --ranges 1:100 101:200 201:300 301:400 401:500
done
```
Extract Tri5 genes promotors regions

```bash

cat gene_pred/codingquarry/F.venenatum/WT_minion/final/Tri5_genes.txt | sed "s/.t.//g" > analysis/meme/promotor_regions/F.venenatum/WT_minion/Tri5_headers.txt

OutDir=analysis/meme/promotor_regions/F.venenatum/WT_minion/tri
mkdir -p $OutDir

for Upstream in $(ls analysis/meme/promotor_regions/F.venenatum/WT_minion/*.upstream*.fasta); do
Region=$(basename ${Upstream%.fasta} | sed 's/promotor_regions.upstream//g')
mkdir $OutDir/$Region
RegionPromotors=$OutDir/$Region/F.venenatum_WT_${Region}_promotors.fa
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
$ProgDir/extract_from_fasta.py --fasta $Upstream --headers analysis/meme/promotor_regions/F.venenatum/WT_minion/Tri5_headers.txt > $RegionPromotors
done
```

```bash
srun --partition medium --mem 10G --cpus-per-task 10 --pty bash
conda activate meme_v4

for Query in $(ls analysis/meme/promotor_regions/F.venenatum/WT_minion/tri/*/*_promotors.fa); do
Region=$(basename ${Query%.fa} | sed 's/.*upstream//g' | sed "s/_prom.*//g")
echo $Region
OutDir=$(dirname $Query)
mkdir -p $OutDir/meme
meme $Query -dna -oc $OutDir/meme -nostatus -mod zoops -nmotifs 5 -minw 6 -maxw 20 -revcomp
ls $OutDir/meme/meme.txt
mast $OutDir/meme/meme.xml $Query -oc $OutDir/meme -nostatus
mv $OutDir/meme/mast.txt $OutDir/meme/${Region}_mast.txt
mv $OutDir/meme/mast.html $OutDir/meme/${Region}_mast.html
done

for Query in $(ls analysis/meme/promotor_regions/F.venenatum/WT_minion/tri/*/*_promotors.fa); do
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