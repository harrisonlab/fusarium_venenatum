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
for File in $(ls gene_pred/braker/*/*_UTR/*/augustus.hints.gff3); do
	Strain=$(echo $File | rev | cut -d '/' -f3 | rev | sed 's/_UTR//g')
	Organism=$(echo $File | rev | cut -d '/' -f4 | rev)
	echo "$Organism - $Strain"
	echo "number of genes:"
	cat $File | grep -v '#' | grep -w 'gene' | wc -l
	echo "number of genes with predicted UTRs"
	cat ${File%.gff3}_utr.gff | grep -v '#' | grep -w 'gene' | wc -l

	getAnnoFasta.pl $File
	OutDir=$(dirname $File)
	echo "##gff-version 3" > $OutDir/augustus_extracted.gff
	cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
done
```

```
  V.dahliae - 12008
  number of genes:
  9602
  number of genes with predicted UTRs
  9196
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
  Fv_genes.gbk

```

## Stand alone run of CASSIS


```bash
screen -a
ssh compute02
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

for GeneID in $(cat anchor_genes.txt | head -n2 | tail -n1); do
  echo $GeneID
  mkdir out/$GeneID
  cassis \
    --annotation cassis.tsv \
    --genome $Assembly \
    --anchor $GeneID \
    --dir out/$GeneID \
    --mismatches 0 \
    -v \
    --prediction \
    --num-cpus 40
done
```

```
Cannot read from binding sites file "out/g717.t1/g717.t1/fimo/+0_-3/fimo.txt".
No such file or directory at script/cassis.pl line 1519.
  +0_-3   (base) armita@compute02:~/tmp/cassis_Fv$ ls out/g717.t1/g717.t1/fimo/+0_-3/fimo.txt
ls: cannot access 'out/g717.t1/g717.t1/fimo/+0_-3/fimo.txt': No such file or directory
```
