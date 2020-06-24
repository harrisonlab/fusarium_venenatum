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



```bash
for Assembly in $(ls rassembly/previous_versions/F.venenatum/WT/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev) 
Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/braker_v2/$Organism/$Strain/
AcceptedHits=alignment/star/F.venenatum/WT/concatenated/concatenated_258647.bam
GeneModelName="$Organism"_"$Strain"_braker
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
sbatch $ProgDir/braker_fungi_v2.sh $Assembly $OutDir $AcceptedHits $GeneModelName
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