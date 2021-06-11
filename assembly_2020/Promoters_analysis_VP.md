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



Extract promotor regions:
Promotor regions upstream of all genes were extracted:

```bash
# This will extract promoter regions from the start codon of a gene. If two genes are close together, this promoter might include exons from the other gene.
for GeneGff in $(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.gff3); do
Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
Assembly=$(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
OutDir=analysis/promoters_VP/promoters_regions/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Promoter_analysis
$ProgDir/extract_promotor.pl --fasta $Assembly --gff $GeneGff --prefix $OutDir/${Strain} --ranges 1:1000 1:2000 1:3000
done
```


```bash
# Create fasta files of RxLR upstream regions
for Upstream in $(ls analysis/promoters_VP/promoters_regions/F.venenatum/WT_minion/*.upstream*000.fasta); do
Region=$(basename ${Upstream%.fasta} | sed 's/promotor_regions.upstream//g')
mkdir $OutDir/$Region
RegionPromotors=$OutDir/$Region/F.venenatum_WT_${Region}_promotors.fa
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
$ProgDir/extract_from_fasta.py --fasta $Upstream --headers $OutDir/tri_headers.txt > $RegionPromotors
done
```

```bash

for Query in $(ls analysis/promoters_VP/promoters_regions/F.venenatum/WT_minion/*000/*_promotors.fa); do
Region=$(basename ${Query%.fa} | sed 's/.*upstream//g' | sed "s/_prom.*//g")
echo $Region
OutDir=$(dirname $Query)
mkdir -p $OutDir/meme
meme $Query -dna -oc $OutDir/meme -nostatus -mod anr -nmotifs 20 -minw 6 -maxw 20 -evt 1.0e+005 -objfun classic -revcomp
ls $OutDir/meme/meme.txt
mast $OutDir/meme/meme.xml $Query -oc $OutDir/meme -nostatus
mv $OutDir/meme/mast.txt $OutDir/meme/${Region}_mast.txt
mv $OutDir/meme/mast.html $OutDir/meme/${Region}_mast.html
done

# Same parameters used for the promoters from cassis
for Query in $(ls analysis/promoters_VP/promoters_regions/F.venenatum/WT_minion/*000/*_promotors.fa); do
Region=$(basename ${Query%.fa} | sed 's/.*upstream//g' | sed "s/_prom.*//g")
echo $Region
OutDir=$(dirname $Query)
mkdir -p $OutDir/meme_maxw12
meme $Query -dna -oc $OutDir/meme_maxw12 -nostatus -mod anr -nmotifs 20 -minw 6 -maxw 12 -evt 1.0e+005 -revcomp
ls $OutDir/meme_maxw12/meme.txt
mast $OutDir/meme_maxw12/meme.xml $Query -oc $OutDir/meme_maxw12 -nostatus
mv $OutDir/meme_maxw12/mast.txt $OutDir/meme_maxw12/${Region}_mast.txt
mv $OutDir/meme_maxw12/mast.html $OutDir/meme_maxw12/${Region}_mast.html
done

```
```bash


UpstreamFa=$(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.upstream3000.fasta)
OutDir=analysis/promoters_VP/promoters_regions/F.venenatum/WT_minion/tri_upstream3000
mkdir -p $OutDir
cat $UpstreamFa | sed -n "/g6426/,/g6438/p;/g6438/q" | head -n-1 > $OutDir/tri_cluster_upstream3000.fasta

for Query in $(ls analysis/promoters_VP/promoters_regions/F.venenatum/WT_minion/tri_upstream3000/tri_cluster_upstream3000.fasta); do
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


  #Create a random sample of a 100 promoter control sequences

  scripts=/home/sobczm/bin/popgen/clock/motif_discovery
  # qsub $scripts/sub_fasta_subsample.sh gene_pred/final/F.venenatum/WT_ncbi/final/final_genes_combined.upstream3000.fasta 100
  qsub $scripts/sub_fasta_subsample.sh gene_pred/final/F.venenatum/WT_ncbi/final/final_genes_combined.upstream3000.fasta 13

  OutDir=analysis/meme/promotor_regions/F.venenatum/WT/tri
  # qsub $scripts/sub_ame.sh $OutDir/tri_cluster_upstream3000.fasta final_genes_combined.upstream3000_random_100.fasta min_ace YNAGGCC
  qsub $scripts/sub_ame.sh $OutDir/tri_cluster_upstream3000.fasta final_genes_combined.upstream3000_random_13.fasta Zn_finger GTGA
  mv *_vs_Zn_finger $OutDir/.




cat final_genes_appended_renamed.upstream3000.fasta | sed 's/_upstream3000//g' > final_genes_appended_renamed.upstream3000_edited.fasta
# Create fasta files of RxLR upstream regions
for Upstream in $(ls *edited.fasta); do
mkdir upstream3000
RegionPromotors=upstream3000/F.venenatum_WT_upstream3000_promotors.fa
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
$ProgDir/extract_from_fasta.py --fasta $Upstream --headers tri_headers.txt > $RegionPromotors
done
```

# Fimo in upstream 3000
```bash
for Query in $(ls upstream3000/F.venenatum_WT_upstream3000_promotors.fa); do
iupac2meme TNAGGCCT > motif1.txt
mkdir -p upstream3000/motif1
fimo -thresh 0.00006 -oc upstream/motif1 motif1.txt $Query
done

for Query in $(ls upstream3000/F.venenatum_WT_upstream3000_promotors.fa); do
iupac2meme YNAGCCC > motif2.txt
fimo -thresh 0.0006 -oc upstream3000/motif2 motif2.txt $Query
done

for Query in $(ls upstream3000/F.venenatum_WT_upstream3000_promotors.fa); do
iupac2meme GTGA > motif3.txt
fimo -thresh 0.006 -oc upstream3000/motif3 motif3.txt $Query
done

for Query in $(ls upstream3000/F.venenatum_WT_upstream3000_promotors.fa); do
iupac2meme SYGGRG > motif4.txt
fimo -thresh 0.0006 -oc upstream3000/motif4 motif4.txt $Query
done
```

# Fimo in cassis promoters
```bash
for Query in $(ls cassis_promoters/tri_promoters.fa); do
#iupac2meme TNAGGCCT > motif1.txt
#mkdir -p cassis_promoters/motif1
fimo -thresh 0.0006 -oc cassis_promoters/motif1 motif1.txt $Query
done

for Query in $(ls cassis_promoters/tri_promoters.fa); do
#iupac2meme YNAGCCC > motif2.txt
fimo -thresh 0.006 -oc cassis_promoters/motif2 motif2.txt $Query
done

for Query in $(ls cassis_promoters/tri_promoters.fa); do
#iupac2meme GTGA > motif3.txt
fimo -thresh 0.006 -oc cassis_promoters/motif3 motif3.txt $Query
done

for Query in $(ls cassis_promoters/tri_promoters.fa); do
#iupac2meme SYGGRG > motif4.txt
fimo -thresh 0.006 -oc cassis_promoters/motif4 motif4.txt $Query
done
```
