# F.gramminearium

Informatics on the downloaded genome of Fusarium gramminearium

the assembly was downloaded from:
http://fungi.ensembl.org/Fusarium_graminearum/Info/Index?db=core
using assembly version RR1

Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/fusarium_venenatum

```bash
  cd /home/groups/harrisonlab/project_files/fusarium_venenatum
```

# Data download

Assemblies and gene models were downloaded for Fg from Emsembl at
http://fungi.ensembl.org/Fusarium_graminearum/Info/Index

these data were copied onto the cluster in the following locations:

```bash
  Organism=F.graminearum
  Strain=PH1
  mkdir -p assembly/external_group/$Organism/$Strain
  mkdir -p assembly/external_group/$Organism/$Strain/gff
```

Data was unzipped prior to further analysis:

```bash
  cd assembly/external_group/F.graminearum/PH1
  for Dir in $(ls -d *); do
    cd $Dir;
    gunzip *.gz;
    cd ../;
  done
  cd ../../../../
```

# Collecting assembly statistics


```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/external_group/F.graminearum/PH1/dna/Fusarium_graminearum.RR1.dna.toplevel.fa); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    echo "$Organism - $Strain"
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```


## Assessing the Gene space in the assembled genome:

```bash
for Assembly in $(ls assembly/external_group/F.graminearum/PH1/dna/Fusarium_graminearum.RR1.dna.toplevel.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
# BuscoDB="Fungal"
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```

```bash
	for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt); do  
		echo $File;
		cat $File | grep -e '(C)' -e 'Total';
	done
```
<!--
# identifying the % of the assembly that has been repeatmasked:

```bash

Unmasked=$(ls assembly/external_group/F.graminearum/PH1/dna/Fusarium_graminearum.RR1.dna.toplevel.fa)
cat $Unmasked | grep -v '>' | grep -o 'N' | wc -l
Hardmasked=$(ls assembly/external_group/F.graminearum/PH1/dna/Fusarium_graminearum.RR1.dna_rm.toplevel.fa)
cat $Hardmasked | grep -v '>' | grep -o 'N' | wc -l
``` -->


# Repeatmasking

Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The best assembly was used to perform repeatmasking

```bash
  for Assembly in $(ls assembly/external_group/F.graminearum/PH1/dna/Fusarium_graminearum.RR1.dna.toplevel.fa); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
    BestAss=assembly/spades/F.venenatum/WT/ncbi_edits/contigs_min_500bp_renamed.fasta
    OutDir=repeat_masked/F.venenatum/WT/illumina_assembly_ncbi
    qsub $ProgDir/rep_modeling.sh $BestAss $OutDir
    qsub $ProgDir/transposonPSI.sh $BestAss $OutDir
  done
```  

<!--
** % bases maked by repeatmasker: 4.75%**

** % bases masked by transposon psi: 4.19% **

The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and hardmasked files.

```bash

for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep -w 'WT' | grep 'ncbi'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa | grep -w 'WT' | grep 'ncbi'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```

```bash
for RepDir in $(ls -d repeat_masked/F.*/*/* | grep -w 'WT' | grep 'ncbi'); do
Strain=$(echo $RepDir | rev | cut -f2 -d '/' | rev)
Organism=$(echo $RepDir | rev | cut -f3 -d '/' | rev)  
RepMaskGff=$(ls $RepDir/*_contigs_hardmasked.gff)
TransPSIGff=$(ls $RepDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
printf "$Organism\t$Strain\n"
# printf "The number of bases masked by RepeatMasker:\t"
sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
# printf "The number of bases masked by TransposonPSI:\t"
sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
# printf "The total number of masked bases are:\t"
cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
echo
done
```

F.venenatum	WT
302604
144657
438768 -->


## Assessing the Gene space in predicted transcriptomes:

```bash
for Assembly in $(ls assembly/external_group/F.graminearum/PH1/cdna/Fusarium_graminearum.RR1.cdna.all.fa); do
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
	for File in $(ls gene_pred/busco/*/*/genes/*/short_summary_*.txt); do  
		echo $File;
		cat $File | grep -e '(C)' -e 'Total';
	done
```


## Identification of SecMet clusters in Fv from Fg

Using Secmet clusters and gene Ids as described in Sieber 2014.


This required download of the Fg 2014 genome from ncbi. This is the same assembly / gene models as retrieved from ensembl but with different geneIDs


```bash
mkdir -p assembly/external_group/F.graminearum/PH1/ncbi
mkdir -p analysis/secondary_metabolites/homologs/F.gramminearum/PH1
```

Data was downlaoded to my local machine and uploaded to the cluster.

```bash
# King V1 PH1 assembly (2014 assembly with gene models updated in 2017)
scp /Volumes/GGB/Harrison/Projects/AHDB\ Fusarium\ C-300009/Science/WP1.3\ Development\ of\ whole\ AMPSeq\ for\ Fus\ communities/2019/Fo\ Lactucae/Fg\ genome/PH1\ v1/GCA_900044135.1_GZPH1RResV1_genomic.* cluster:/home/groups/harrisonlab/project_files/fusarium_venenatum/assembly/external_group/F.graminearum/PH1/ncbi/.
# Broad V3 PH1 assemebly (2003 assembly with gene models updated in 2014)
scp /Volumes/GGB/Harrison/Projects/BBSRC\ IPA\ QUORN\ C300045/Science/Obj1\ Reference\ generation/2019/Presence\ of\ Fg\ secondary\ metabolite\ clusters/Fg\ genome/PH1\ v3/GCF_000240135.3_ASM24013v3_genomic* cluster:/home/groups/harrisonlab/project_files/fusarium_venenatum/assembly/external_group/F.graminearum/PH1/ncbi/.

```
A file of Sieber 2014 gene clusters was uploaded to the cluster:

```bash
scp /Volumes/GGB/Harrison/Projects/AHDB\ Fusarium\ C-300009/Science/WP1.3\ Development\ of\ whole\ AMPSeq\ for\ Fus\ communities/2019/Fo\ Lactucae/Fg\ genome/sieber_clusters.tsv cluster:/home/groups/harrisonlab/project_files/fusarium_venenatum/analysis/secondary_metabolites/homologs/F.gramminearum/PH1/.
```

```bash
cd /home/groups/harrisonlab/project_files/fusarium_venenatum

sed -i "s/ - /\t/g" analysis/secondary_metabolites/homologs/F.gramminearum/PH1/sieber_clusters.tsv

# For King v1 genes
Gff=$(ls assembly/external_group/F.graminearum/PH1/ncbi*/*.gff | grep 'king_v1')
FgClusters=$(ls analysis/secondary_metabolites/homologs/F.gramminearum/PH1/sieber_clusters.tsv)
OutDir=$(dirname $FgClusters)/king_v1
mkdir $OutDir

rm $OutDir/Fg_cluster_locations.tsv
rm $OutDir/Fg_cluster_locations.gff

for Cluster in $(cat $FgClusters | tail -n +2 | cut -f1); do
GeneStart=$(cat $FgClusters | grep -w "$Cluster" | cut -f2)
GeneStop=$(cat $FgClusters | grep -w "$Cluster" | cut -f3)
# printf "$Cluster\t$GeneStart\t$GeneStop\n"
LocationStart=$(cat $Gff | grep -w $GeneStart | head -n1 | cut -f1,4)
LocationStop=$(cat $Gff | grep -w $GeneStop | head -n1 | cut -f1,5)
printf "$Cluster\t$LocationStart\t$LocationStop\n" >> $OutDir/Fg_cluster_locations.tsv
Contig=$(echo $LocationStart | cut -f1 -d ' ')
Start=$(echo $LocationStart | cut -f2 -d ' ')
Stop=$(echo $LocationStop | cut -f2 -d ' ')
if [ "$Start" ] && [ $Stop ]; then
  if [[ $Stop -gt $Start ]]; then
    printf "$Contig\tNIAB_EMR\tFg_cluster\t$Start\t$Stop\t.\t.\t.\tID=$Cluster;\n" >> $OutDir/Fg_cluster_locations.gff
  else
    printf "$Contig\tNIAB_EMR\tFg_cluster\t$Stop\t$Start\t.\t.\t.\tID=$Cluster;\n" >> $OutDir/Fg_cluster_locations.gff
  fi
fi
done


# For Broad v3 genes
Gff=$(ls assembly/external_group/F.graminearum/PH1/ncbi*/*.gff* | grep 'broad_v3')
FgClusters=$(ls analysis/secondary_metabolites/homologs/F.gramminearum/PH1/sieber_clusters.tsv)
OutDir=$(dirname $FgClusters)/broad_v3
mkdir $OutDir

rm $OutDir/Fg_cluster_locations.tsv
rm $OutDir/Fg_cluster_locations.gff

for Cluster in $(cat $FgClusters | tail -n +2 | cut -f1); do
GeneStart=$(cat $FgClusters | grep -w "$Cluster" | cut -f2)
GeneStop=$(cat $FgClusters | grep -w "$Cluster" | cut -f3)
# printf "$Cluster\t$GeneStart\t$GeneStop\n"
LocationStart=$(cat $Gff | grep -w $GeneStart | head -n1 | cut -f1,4)
LocationStop=$(cat $Gff | grep -w $GeneStop | head -n1 | cut -f1,5)
printf "$Cluster\t$LocationStart\t$LocationStop\n" >> $OutDir/Fg_cluster_locations.tsv
Contig=$(echo $LocationStart | cut -f1 -d ' ')
Start=$(echo $LocationStart | cut -f2 -d ' ')
Stop=$(echo $LocationStop | cut -f2 -d ' ')
if [ "$Start" ] && [ $Stop ]; then
  if [[ $Stop -gt $Start ]]; then
    printf "$Contig\tNIAB_EMR\tFg_cluster\t$Start\t$Stop\t.\t.\t.\tID=$Cluster;\n" >> $OutDir/Fg_cluster_locations.gff
  else
    printf "$Contig\tNIAB_EMR\tFg_cluster\t$Stop\t$Start\t.\t.\t.\tID=$Cluster;\n" >> $OutDir/Fg_cluster_locations.gff
  fi
fi
done
```

Some clusters were missed in each analysis:
```bash
cat analysis/secondary_metabolites/homologs/F.gramminearum/PH1/king_v1/Fg_cluster_locations.tsv | grep -v -P "\t.+\t.+\t.+\t" | cut -f1 | wc -l
cat analysis/secondary_metabolites/homologs/F.gramminearum/PH1/broad_v3/Fg_cluster_locations.tsv | grep -v -P "\t.+\t.+\t.+\t" | cut -f1 | wc -l

cat analysis/secondary_metabolites/homologs/F.gramminearum/PH1/*/Fg_cluster_locations.tsv | grep -v -P "\t.+\t.+\t.+\t" | cut -f1 | sort | uniq -d
```

```
C04
C12
C15
C17
C33
C35
C39
C42
C46
C48
C53
C54
C61
C64
C65
```

To associate Fg gene clusters with Fv, Fv proteins were blasted
against Fg genomes and where they intersected with Fg gene clusters
then


```bash
for Assembly in $(ls assembly/external_group/F.graminearum/PH1/ncbi*/*.fna | grep 'broad'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Version=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
echo $Assembly
Query=gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.gene.fasta
OutDir=analysis/secondary_metabolites/$Organism/$Strain/$Version
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query dna $Assembly $OutDir
done
```

<!-- Manual searches

```bash
qlogin -pe smp 14
cd /home/groups/harrisonlab/project_files/fusarium_venenatum
for Assembly in $(ls assembly/external_group/F.graminearum/PH1/ncbi*/*.fna); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Version=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
echo $Assembly
Query=$(ls gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.gene.fasta)
OutDir=analysis/FTF/$Organism/$Strain/${Version}_2
mkdir -p $OutDir
dbFasta=$(ls /home/groups/harrisonlab/phibase/v4.4/phi_accessions.fa)
dbType="nucl"
Prefix="${Strain}_${Version}"
Eval="1e-30"
makeblastdb -in $Assembly -input_type fasta -dbtype $dbType -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
tblastx -num_threads 14 -db $OutDir/$Prefix.db -query $Query -outfmt 6 -num_alignments 1 -out $OutDir/${Prefix}_hits.txt -evalue $Eval
cat $OutDir/${Prefix}_hits.txt | cut -f1,2 | sort | uniq > $OutDir/${Prefix}_hits_headers.txt
done > log.txt
``` -->

BLAST hits were converted to Gff annotations and intersected with gene models:

```bash
for BlastHits in $(ls analysis/secondary_metabolites/F.graminearum/PH1/*/PH1_final_genes_appended_renamed.gene.fasta_hits.csv); do
Strain=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
Organism=$(echo $BlastHits | rev | cut -f4 -d '/' | rev)
Version=$(echo $BlastHits | rev | cut -f2 -d '/' | rev | sed 's/ncbi_//g')
# OutDir=analysis/FTF/$Organism/$Strain
HitsGff=$(echo $BlastHits | sed  's/.csv/.gff/g')
OutDir=$(dirname $BlastHits)
Column2=Fv_gene_homolog
NumHits=1
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
# $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff


GffClusters=$(ls analysis/secondary_metabolites/homologs/F.gramminearum/PH1/$Version/Fg_cluster_locations.gff)
bedtools intersect -wo -a $HitsGff -b $GffClusters > $OutDir/"$Strain"_SecMet_hits_intersected.bed
echo "$Version"
for Cluster in $(cat $OutDir/"$Strain"_SecMet_hits_intersected.bed | cut -f18 | sed 's/ID=//g' | tr -d ';' | sort | uniq); do
  mkdir -p $OutDir/$Cluster
  cat $OutDir/"$Strain"_SecMet_hits_intersected.bed | grep "$Cluster" | cut -f1-9 | sed "s/Fv_gene_homolog/${Cluster}_homolog/g" > $OutDir/$Cluster/${Strain}_SecMet_${Cluster}_hits_intersected.bed
  cat $OutDir/$Cluster/${Strain}_SecMet_${Cluster}_hits_intersected.bed | cut -f9 | sed "s/\"ID=//g" | sed "s/_Blast.*/.t1/g" > $OutDir/$Cluster/${Strain}_SecMet_${Cluster}_hits_intersected_headers.txt
  Gff=$(ls gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.gff3)
  ProgDir=/home/passet/git_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_gff_for_sigP_hits.pl $OutDir/$Cluster/${Strain}_SecMet_${Cluster}_hits_intersected_headers.txt $Gff ${Cluster}_homolog ID > $OutDir/$Cluster/${Strain}_${Version}_${Cluster}_hits.gff
  printf "$Cluster\t"
  cat $OutDir/$Cluster/${Strain}_SecMet_${Cluster}_hits_intersected_headers.txt | sed -r "s/$/,/g" | tr -d "\n"
echo ""
done
done


for GffClusters in $(ls analysis/secondary_metabolites/homologs/F.gramminearum/PH1/*/Fg_cluster_locations.gff); do
  echo $GffClusters
  cat $GffClusters | awk '{print $9, $1, $4, $5}' | sed 's/ID=//g' | tr -d ";"
done
```
