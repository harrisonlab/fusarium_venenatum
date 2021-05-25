Extract promotor regions:
Promotor regions upstream of all genes were extracted:

```bash
for GeneGff in $(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.gff3); do
Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
Assembly=$(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
OutDir=analysis/promoters_VP/promoters_regions/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Promoter_analysis
$ProgDir/extract_promotor.pl --fasta $Assembly --gff $GeneGff --prefix $OutDir/${Strain} --ranges 1:100 101:200 201:300 301:400 401:500
done
```
## Extract cluster intergenic regions

These Tri5 genes were extracted:

```bash

OutDir=analysis/meme/promotor_regions/F.venenatum/WT/tri
mkdir -p $OutDir
TriGenes="g6426 g6427 g6428 g6429 g6430 g6431 g6432 g6433 g6434 g6435 g6436 g6437"
printf "$TriGenes" | sed "s/ /\n/g" | sed "s/.t.//g" > $OutDir/tri_headers.txt


# Create fasta files of RxLR upstream regions
for Upstream in $(ls analysis/promoters_VP/promoters_regions/F.venenatum/WT_minion/*.upstream*.fasta); do
Region=$(basename ${Upstream%.fasta} | sed 's/promotor_regions.upstream//g')
mkdir $OutDir/$Region
RegionPromotors=$OutDir/$Region/F.venenatum_WT_${Region}_promotors.fa
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
$ProgDir/extract_from_fasta.py --fasta $Upstream --headers $OutDir/tri_headers.txt > $RegionPromotors
done
```
```bash

for Query in $(ls analysis/promoters_VP/promoters_regions/F.venenatum/WT_minion/*100/*_promotors.fa); do
Region=$(basename ${Query%.fa} | sed 's/.*upstream//g' | sed "s/_prom.*//g")
echo $Region
OutDir=$(dirname $Query)
mkdir -p $OutDir/meme
meme $Query -dna -oc $OutDir/meme -nostatus -mod anr -nmotifs 5 -minw 6 -maxw 20 -objfun classic -revcomp
ls $OutDir/meme/meme.txt
mast $OutDir/meme/meme.xml $Query -oc $OutDir/meme -nostatus
mv $OutDir/meme/mast.txt $OutDir/meme/${Region}_mast.txt
mv $OutDir/meme/mast.html $OutDir/meme/${Region}_mast.html
done

```

No significant motif founded. I need larger regions