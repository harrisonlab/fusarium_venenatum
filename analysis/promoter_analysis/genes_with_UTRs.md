# Genes with UTRs

Gene models were repredicted using an updated version of Braker
on the new cluster that allows UTR prediction.

This location is in the /data2 area of the old cluster:

```bash
  ProjDir=/oldhpc/data/scratch/armita/fusarium_venenatum
  cd $ProjectDir
```

Gene prediction was run from a screen session with a ssh connection to a worker node:

```bash
screen -a
ssh compute01
cd /oldhpc/data/scratch/armita/fusarium_venenatum


WorkDir=$HOME/tmp/braker_Fv
OldProjDir=/oldhpc/home/groups/harrisonlab/project_files/fusarium_venenatum

Assembly=$(ls $OldProjDir/repeat_masked/F.venenatum/WT/illumina_assembly_ncbi/WT_contigs_softmasked_repeatmasker_TPSI_appended.fa)
Organism="F.venenatum"
Strain="WT"

OutDir=/oldhpc/data/scratch/armita/fusarium_venenatum/fusarium_venenatum/gene_pred/braker/F.venenatum/WT_braker_UTR
AcceptedHits=$(ls $OldProjDir/alignment/star/F.venenatum/WT/concatenated/concatenated.bam)
GeneModelName="$Organism"_"$Strain"_braker
CurDir=$PWD

mkdir -p $WorkDir
cd $WorkDir

cp $Assembly assembly.fa
cp $AcceptedHits alignedRNA.bam

braker.pl \
  --cores 40 \
  --overwrite \
  --fungus \
  --UTR=on \
  --gff3 \
  --softmasking on \
  --species=$GeneModelName \
  --genome="assembly.fa" \
  --bam="alignedRNA.bam"

mkdir -p $CurDir/$OutDir
cp -r braker/* $CurDir/$OutDir/.

rm -r $WorkDir

```
