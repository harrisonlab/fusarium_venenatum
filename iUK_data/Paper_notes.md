# Validation of results for Paper 1

Total number of secondary metabolite clusters:

```bash
cd WT_minion
for Assembly in $(ls assembly/repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/medaka_contigs_unmasked.fa); do
Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
# OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
# mkdir -p $OutDir
GeneGff=$(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.gff3)
AntismashClusters=$(ls analysis/secondary_metabolites/antismash/WT_minion_VP/WT_antismash_results_secmet_clusters.gff)
SmurfClusters=$(ls analysis/secondary_metabolites/smurf/F.venenatum/WT_minion/*lusters.gff)
echo "Total number of Antismash clusters"
cat $AntismashClusters | wc -l
echo "Total number of SMURF clusters"
cat $SmurfClusters | wc -l
echo "number of Antismash clusters intersecting Smurf clusters"
bedtools intersect -a $AntismashClusters -b $SmurfClusters | wc -l
echo "number of Antismash clusters not intersecting Smurf clusters"
bedtools intersect -v -a $AntismashClusters -b $SmurfClusters | wc -l
echo "number of smurf clusters intersecting antismash clusters"
bedtools intersect -a $SmurfClusters -b $AntismashClusters | wc -l
echo "number of smurf clusters not intersecting antismash clusters"
bedtools intersect -v -a $SmurfClusters -b $AntismashClusters | wc -l
done
```
```
Total number of Antismash clusters
48
Total number of SMURF clusters
23
number of Antismash clusters intersecting Smurf clusters
22
number of Antismash clusters not intersecting Smurf clusters
28
number of smurf clusters intersecting antismash clusters
22
number of smurf clusters not intersecting antismash clusters
9
```

```bash
OutDir=projects/niab/agomez/fusarium_venenatum/WT/analysis/secondary_metabolites/smurf/WT
Prefix="WT"
GeneGff=projects/niab/agomez/fusarium_venenatum/WT/gene_pred/codingquarry/F.venenatum/WT/final/final_genes_appended_renamed.gff3
SmurfClusters=$OutDir/Secondary-Metabolite-Clusters.txt
SmurfBackbone=$OutDir/Backbone-genes.txt
tools/smurf2gff.py --smurf_clusters $SmurfClusters --smurf_backbone $SmurfBackbone > $OutDir/Smurf_clusters.gff
bedtools intersect -wo -a $GeneGff -b $OutDir/Smurf_clusters.gff | grep 'mRNA' | cut -f9,10,12,18 | sed "s/ID=//g" | perl -p -i -e "s/;Parent=g\w+//g" | perl -p -i -e "s/;Notes=.*//g" > $OutDir/"$Prefix"_smurf_secmet_genes.tsv

cd WT
for Assembly in $(ls assembly/WT/WT_contigs_unmasked.fa); do
GeneGff=$(ls gene_pred/codingquarry/F.venenatum/WT/final/final_genes_appended_renamed.gff3)
AntismashClusters=$(ls analysis/secondary_metabolites/antismash/WT/WT_antismash_results_secmet_clusters.gff)
SmurfClusters=$(ls analysis/secondary_metabolites/smurf/WT/*lusters.gff)
echo "Total number of Antismash clusters"
cat $AntismashClusters | wc -l
echo "Total number of SMURF clusters"
cat $SmurfClusters | wc -l
echo "number of Antismash clusters intersecting Smurf clusters"
bedtools intersect -a $AntismashClusters -b $SmurfClusters | wc -l
echo "number of Antismash clusters not intersecting Smurf clusters"
bedtools intersect -v -a $AntismashClusters -b $SmurfClusters | wc -l
echo "number of smurf clusters intersecting antismash clusters"
bedtools intersect -a $SmurfClusters -b $AntismashClusters | wc -l
echo "number of smurf clusters not intersecting antismash clusters"
bedtools intersect -v -a $SmurfClusters -b $AntismashClusters | wc -l
done
```
```
Total number of Antismash clusters
50
Total number of SMURF clusters
22
number of Antismash clusters intersecting Smurf clusters
5
number of Antismash clusters not intersecting Smurf clusters
45
number of smurf clusters intersecting antismash clusters
5
number of smurf clusters not intersecting antismash clusters
18
```

