# AreA genes analysis

## Extract cluster intergenic regions


```bash

OutDir=analysis/meme/promotor_regions/F.venenatum/WT/AreA
mkdir -p $OutDir
AreAGenes="g1528 g5192 g13141"
printf "$AreAGenes" | sed "s/ /\n/g" | sed "s/.t.//g" > $OutDir/AreA_headers.txt


# Create fasta files of RxLR upstream regions
for Upstream in $(ls analysis/promoters_VP/promoters_regions/F.venenatum/WT_minion/*.upstream*.fasta); do
Region=$(basename ${Upstream%.fasta} | sed 's/promotor_regions.upstream//g')
mkdir $OutDir/$Region
RegionPromotors=$OutDir/$Region/F.venenatum_WT_${Region}_promotors.fa
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
$ProgDir/extract_from_fasta.py --fasta $Upstream --headers $OutDir/AreA_headers.txt > $RegionPromotors
done
```
```bash
for Query in $(ls analysis/meme/promotor_regions/F.venenatum/WT/AreA/*/*_promotors.fa); do
Region=$(basename ${Query%.fa} | sed 's/.*upstream//g' | sed "s/_prom.*//g")
echo $Region
OutDir=analysis/meme/promotor_regions/F.venenatum/WT/AreA/"$Region"
mkdir -p $OutDir/meme
#meme $Query -dna -oc $OutDir/meme -nostatus -mod anr -nmotifs 5 -minw 6 -maxw 20 -objfun classic -revcomp
streme --p $Query -dna -oc $OutDir/streme --objfun de -nmotifs 5 -minw 6 -maxw 20
ls $OutDir/meme/meme.txt
mast $OutDir/meme/meme.xml $Query -oc $OutDir/meme -nostatus
mv $OutDir/meme/mast.txt $OutDir/meme/${Region}_mast.txt
mv $OutDir/meme/mast.html $OutDir/meme/${Region}_mast.html
done
