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

OutDir=$OldProjDir/gene_pred/braker/$Organism/${Strain}_UTR
AcceptedHits=$(ls $OldProjDir/alignment/star/$Organism/$Strain/concatenated/concatenated.bam)
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
  # --stranded=. \
  # --AUGUSTUS_ab_initio \
  --gff3 \
  --softmasking on \
  --species=$GeneModelName \
  --genome="assembly.fa" \
  --bam="alignedRNA.bam"

mkdir -p $OutDir
cp -r braker/* $OutDir/.

rm -r $WorkDir

```

Fasta and gff files were extracted from Braker1 output.

```bash
# for File in $(ls gene_pred/braker/*/*_UTR/*/augustus.hints.gff3); do
for File in $(ls gene_pred/braker/*/*_UTR/*/augustus.hints_utr.gff3); do
	Strain=$(echo $File | rev | cut -d '/' -f3 | rev | sed 's/_UTR//g')
	Organism=$(echo $File | rev | cut -d '/' -f4 | rev)
	echo "$Organism - $Strain"
	echo "number of genes:"
	cat $File | grep -v '#' | grep -w 'gene' | wc -l
	# echo "number of genes with predicted UTRs"
	# cat ${File%.gff3}_utr.gff | grep -v '#' | grep -w 'gene' | wc -l

	getAnnoFasta.pl $File
	OutDir=$(dirname $File)
	echo "##gff-version 3" > $OutDir/augustus_extracted.gff
	cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
done
```

<!-- ```
  V.dahliae - 12008
  number of genes:
  9602
  number of genes with predicted UTRs
  9196
``` -->

```
  F.venenatum - WT
  number of genes:
  11481
```

# Create a conversion table of old to new gene IDs

Run from the old cluster:
```bash
cd /home/groups/harrisonlab/project_files/fusarium_venenatum
NewGff=$(ls gene_pred/braker/*/*_UTR/*/augustus.hints_utr.gff3)
OldGff=$(ls gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.gff3)

OutDir=$(dirname $NewGff)/gene_conversion
mkdir -p $OutDir
cat $OldGff | grep -w 'gene' > $OutDir/old_locations.gff
cat $NewGff | grep -w 'gene' > $OutDir/new_locations.gff

bedtools intersect -loj -s -a $OutDir/old_locations.gff -b $OutDir/new_locations.gff > $OutDir/gene_intersects.gff
```

## Antismash secmet prediction with Cassis

Log into the new cluster:

```bash
screen -a
ssh compute02
WorkDir=~/tmp/antismash_Fv
mkdir $WorkDir
cd $WorkDir
OldProjDir=/oldhpc/home/groups/harrisonlab/project_files/fusarium_venenatum
Assembly=$(ls $OldProjDir/repeat_masked/F.venenatum/WT/illumina_assembly_ncbi/WT_contigs_softmasked_repeatmasker_TPSI_appended.fa)
# Genes=$(ls $OldProjDir/gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.gff3)
Genes=$(ls $OldProjDir/gene_pred/braker/F.venenatum/WT_UTR/*/augustus_extracted.gff)
conda activate emboss
seqret -sequence $Assembly -feature -fformat gff -fopenfile $Genes -osformat genbank -auto -outseq Fv_genes2.gbk
conda deactivate

conda activate antismash
antismash \
  --cpus 40 \
  --taxon fungi \
  --transatpks_da \
  --clusterblast \
  --subclusterblast \
  --knownclusterblast \
  --smcogs \
  --inclusive \
  --cassis \
  --borderpredict \
  --full-hmmer \
  --asf \
  --tta \
  --outputfolder Fv_antismash \
  --verbose \
  Fv_genes2.gbk
```


# Identify promoter regions of genes with UTRs.

Extract region from around UTR -1000 to +50 unless intersecting another gene model

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/promoters
Gff=$(ls gene_pred/braker/F.venenatum/WT_UTR/__braker/augustus.hints_utr.gff3)
Assembly=$(ls repeat_masked/F.venenatum/WT/illumina_assembly_ncbi/WT_contigs_unmasked.fa)
OutDir=analysis/promoters/F.venenatum/WT_UTR
mkdir -p $OutDir
$ProgDir/extract_promoters.py --gff $Gff --fasta $Assembly --prefix $OutDir/WT_promoters
```

```
Average distance of TSS from start codon:
461.675954232
```
