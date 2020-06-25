# Gene prediction MiSeq genome

## Star

Spliced Transcripts Alignment to a Reference. 


```bash
conda install star

for Assembly in $(ls assembly/previous_versions/F.venenatum/WT/*_contigs_unmasked.fa)
do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev) 
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo $Assembly
echo "$Organism - $Strain"
for FileF in $(ls -d ../oldhome/groups/harrisonlab/project_files/quorn/filtered/WTCHG_259732_224.1.fq)
do
FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/1.fq/2.fq/g')
echo $FileF
echo $FileR
Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/.1.fq//g')
OutDir=alignment/star/$Organism/$Strain/$Sample_Name
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
sbatch $ProgDir/star2.sh $Assembly $FileF $FileR $OutDir
done
done
```

Accepted hits .bam file were concatenated and indexed for use for gene model training

```bash
#All alignments. This will be used for Braker
screen -a

BamFiles=$(ls alignment/star/F.venenatum/WT/*/*.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
OutDir=/data/scratch/gomeza/fusarium_venenatum/star/F.venenatum/WT/concatenated
mkdir -p $OutDir
samtools merge -f $OutDir/concatenated.bam $BamFiles

# WT media samples only. This will be used for Codingquarry
BamFiles=$(ls alignment/star/F.venenatum/WT/one_per_media/*/*.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
OutDir=alignment/star/F.venenatum/WT/one_per_media/concatenated
mkdir -p $OutDir
samtools merge -f $OutDir/concatenated.bam $BamFiles
```

## Braker

```bash
for Assembly in $(ls assembly/previous_versions/F.venenatum/WT/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev) 
Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/braker/$Organism/$Strain
AcceptedHits=../../data/scratch/gomeza/fusarium_venenatum/star/F.venenatum/WT/concatenated/concatenated.bam
GeneModelName="$Organism"_"$Strain"_braker
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
sbatch $ProgDir/braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```

## Cufflinks RNA-seq alignments assembler

```bash
    for Assembly in $(ls assembly/previous_versions/F.venenatum/WT/*_contigs_unmasked.fa); do
        Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev) 
        Organism=$(echo $Assembly| rev | cut -d '/' -f3 | rev) 
        echo "$Organism - $Strain"
        OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim
        mkdir -p $OutDir
        AcceptedHits=alignment/star/F.venenatum/WT/one_per_media/concatenated/concatenated.bam
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
        sbatch $ProgDir/cufflinks.sh $AcceptedHits $OutDir
    done
```

## Codingquarry 

```bash
conda activate antismash_py27

for Assembly in $(ls assembly/previous_versions/F.venenatum/WT/*_contigs_unmasked.fa); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev) 
  Organism=$(echo $Assembly| rev | cut -d '/' -f3 | rev) 
  echo "$Organism - $Strain"
  OutDir=gene_pred/codingquarry/$Organism/$Strain
  mkdir -p $OutDir
  GTF=gene_pred/cufflinks/F.venenatum/WT/concatenated_prelim/transcripts.gtf
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
  sbatch $ProgDir/codingquarry2.sh $Assembly $GTF $OutDir
done
```

Additional transcripts predicted by CodingQuarry are added to the final gene models.

```bash
  # The following perl scripts requires the installation of some libraries. Run these commands in a perly environment.
  # Install the required libraries (if any) using cpanm
  # cpanm Bio::Perl

BrakerGff=gene_pred/braker/F.venenatum/WT/augustus.hints.gff3
Strain=$(echo $BrakerGff| rev | cut -d '/' -f2 | rev)
Organism=$(echo $BrakerGff | rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
Assembly=$(ls assembly/previous_versions/F.venenatum/WT/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
CodingQuarryGff=gene_pred/codingquarry/F.venenatum/WT/out/PredictedPass.gff3
PGNGff=gene_pred/codingquarry/F.venenatum/WT/out/PGN_predictedPass.gff3
AddDir=gene_pred/codingquarry/$Organism/$Strain/additional
FinalDir=gene_pred/codingquarry/$Organism/$Strain/final
AddGenesList=$AddDir/additional_genes.txt
AddGenesGff=$AddDir/additional_genes.gff
FinalGff=$AddDir/combined_genes.gff
mkdir -p $AddDir
mkdir -p $FinalDir

  # Create a list with the additional transcripts in CondingQuarry gff (and CQPM) vs Braker gene models
bedtools intersect -v -a $CodingQuarryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
  
  # Creat Gff file with the additional transcripts
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuarryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
	
  # Create a final Gff file with gene features
$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3

# Create fasta files from each gene feature in the CodingQuarry gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary

# Create fasta files from each gene feature in the Braker gff3
cp $BrakerGff $FinalDir/final_genes_Braker.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker

# Combine both fasta files
cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

# Combine both gff3 files
GffBraker=$FinalDir/final_genes_CodingQuary.gff3
GffQuary=$FinalDir/final_genes_Braker.gff3
GffAppended=$FinalDir/final_genes_appended.gff3
cat $GffBraker $GffQuary > $GffAppended

# Check the final number of genes

for DirPath in $(ls -d $FinalDir); do
echo $DirPath;
cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
echo "";
done
```
```
gene_pred/codingquarry/F.venenatum/WT/final
12692
1005
13697
```

  #### Remove duplicate and rename genes.

  ```bash
  for GffAppended in $(ls gene_pred/codingquarry/F.venenatum/WT/final/final_genes_appended.gff3);
  do
    Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    FinalDir=gene_pred/codingquarry/F.venenatum/WT/final
    # Remove duplicated genes
    GffFiltered=gene_pred/codingquarry/F.venenatum/WT/final/filtered_duplicates.gff
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
    $ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
    # Rename genes
    GffRenamed=gene_pred/codingquarry/F.venenatum/WT/final/final_genes_appended_renamed.gff3
    LogFile=gene_pred/codingquarry/F.venenatum/WT/final/final_genes_appended_renamed.log
    $ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
    rm $GffFiltered
    # Create renamed fasta files from each gene feature   
    Assembly=$(ls assembly/previous_versions/F.venenatum/WT/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
    $ProgDir/gff2fasta.pl $Assembly $GffRenamed gene_pred/codingquarry/F.venenatum/WT/final/final_genes_appended_renamed
    # The proteins fasta file contains * instead of Xs for stop codons, these should be changed
    sed -i 's/\*/X/g' gene_pred/codingquarry/F.venenatum/WT/final/final_genes_appended_renamed.pep.fasta
  done 
```

## Busco

```bash
  conda activate BUSCO
  for Fasta in $(ls gene_pred/codingquarry/F.venenatum/WT/final/final_genes_appended_renamed.gene.fasta); do
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
    BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
    OutDir=$(dirname $Fasta)/busco_sordariomycetes_obd10
    sbatch $ProgDir/busco.sh $Fasta $BuscoDB $OutDir
  done
```

## Interproscan

Interproscan was used to give gene models functional annotations.


```bash
# This command will split your gene fasta file and run multiple interproscan jobs.
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  for Genes in $(ls gene_pred/codingquarry/F.venenatum/WT/final/final_genes_appended_renamed.pep.fasta); do
    echo $Genes
    $ProgDir/interproscan.sh $Genes
  done 2>&1 | tee -a interproscan_submisison.log
```







## SwissProt

```bash
for Proteome in $(ls gene_pred/codingquarry/F.venenatum/WT/final/final_genes_appended_renamed.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/swissprot/$Organism/$Strain
SwissDbDir=../dbUniprot/swissprot_2020_June
SwissDbName=uniprot_sprot.db
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
sbatch $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
done
```

```bash
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
CurPath=$PWD
for Proteome in $(ls gene_pred/codingquarry/F.venenatum/WT/final/final_genes_appended_renamed.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
SplitDir=gene_pred/final_genes_split/$Organism/$Strain
mkdir -p $SplitDir
BaseName="$Organism""_$Strain"_final_preds
$ProgDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
for File in $(ls $SplitDir/*_final_preds_*); do
sbatch $ProgDir/pred_signalP.sh $File signalp-4.1
done
done
```