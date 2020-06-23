## Star

Spliced Transcripts Alignment to a Reference. 

### Requirements

```bash
conda install star
```

### Typical run

```bash

for Assembly in $(ls assembly/previous_versions/F.venenatum/WT/*_contigs_unmasked.fa)
do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev) 
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo $Assembly
echo "$Organism - $Strain"
for FileF in $(ls -d ../oldhome/groups/harrisonlab/project_files/quorn/filtered/WTCHG_258647_224.1.fq)
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

BamFiles=$(ls alignment/star/F.venenatum/WT/*/*.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
OutDir=alignment/star/F.venenatum/WT/concatenated
mkdir -p $OutDir
samtools merge -f $OutDir/concatenated_258647.bam $BamFiles

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