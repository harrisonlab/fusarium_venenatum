# Fusarium venenatum summary

The updated files are stored in gomez_WD

## Build annotation tables

```bash
    for GeneGff in $(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.gff3); do
        Strain=WT_minion
        Organism=F.venenatum
        Assembly=$(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
        TFs=$(ls analysis/transcription_factors/F.venenatum/WT_minion/WT_minion_TF_domains.tsv )
        InterPro=$(ls gene_pred/interproscan/F.venenatum/WT_minion/WT_minion_interproscan.tsv)
        Antismash=$(ls analysis/secondary_metabolites/antismash/F.venenatum/WT_minion/WT_antismash_results_secmet_genes_corrected.tsv)
        SwissProt=$(ls gene_pred/swissprot/F.venenatum/WT_minion/swissprot_vJun2020_tophit_parsed.tbl)
        OutDir=analysis/annotation_tables
        mkdir -p $OutDir
        GeneFasta=$(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta)
        Dir1=$(ls -d analysis/RNAseq/salmon/F.venenatum/WT_minion/DeSeq2/Contrasts)
        DEG_Files=$(ls \
        $Dir1/RH2_vs_RH1.txt \
        $Dir1/RH3_vs_RH1.txt \
        $Dir1/RH4_vs_RH1.txt \
        $Dir1/RH5_vs_RH1.txt \
        $Dir1/RH6_vs_RH1.txt \
        $Dir1/RH7_vs_RH1.txt \
        $Dir1/RH8_vs_RH1.txt \
        $Dir1/RH7_vs_RH5.txt \
        | sed -e "s/$/ /g" | tr -d "\n")
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Annotation_tables
        $ProgDir/build_annot_RNAseq.py --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --DEG_files $DEG_Files --TFs $TFs --InterPro $InterPro --Antismash $Antismash --Swissprot $SwissProt > $OutDir/"$Strain"_gene_table.tsv
    done
```

## Investigate enriched functional annotations in DEGs vs all genes

```bash
    OutDir=analysis/enrichment/F.venenatum/WT_minion/whole_genome
    mkdir -p $OutDir
    InterProTSV=gene_pred/interproscan/F.venenatum/WT_minion/WT_minion_interproscan.tsv
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
    $ProgDir/GO_table.py --interpro $InterProTSV > $OutDir/experiment_all_gene_GO_annots.tsv
    cat $OutDir/experiment_all_gene_GO_annots.tsv | sed 's/.t.*//g' > $OutDir/temp1.tsv
    cat $OutDir/experiment_all_gene_GO_annots.tsv | cut -f2 > $OutDir/temp2.tsv
    paste $OutDir/temp1.tsv $OutDir/temp2.tsv > $OutDir/experiment_all_gene_GO_annots_geneid.tsv
    rm $OutDir/temp1.tsv
    rm $OutDir/temp2.tsv
```

```bash
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
AnnotTable=analysis/annotation_tables/WT_minion_gene_table.tsv
DEGs=analysis/RNAseq/salmon/F.venenatum/WT_minion/DeSeq2/Contrasts/RH2_vs_RH1_DEGs.txt
AllGenes=$OutDir/Allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/RH2vRH1all.txt
Set2Genes=$OutDir/RH2vRH1all2.txt
AllGenes=$OutDir/t1_down_genes.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes
$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt
```


 ls analysis/Coexpression/Clust_Results_m70xq10z2qz49uwe/Clusters_Objects.tsv | cut -f1 | sed "1,2d" > Cluster1.txt
 ls analysis/Coexpression/Clust_Results_m70xq10z2qz49uwe/Clusters_Objects.tsv | cut -f2 | sed "1,2d" > Cluster2.txt

```bash
OutDir=analysis/Coexpression/enrichment/F.venenatum/WT_minion/whole_genome/C0
mkdir -p $OutDir
InterProTSV=gene_pred/interproscan/F.venenatum/WT_minion/WT_minion_interproscan.tsv
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
$ProgDir/GO_table.py --interpro $InterProTSV > $OutDir/experiment_all_gene_GO_annots.tsv
cat $OutDir/experiment_all_gene_GO_annots.tsv | sed 's/.t.*//g' > $OutDir/temp1.tsv
cat $OutDir/experiment_all_gene_GO_annots.tsv | cut -f2 > $OutDir/temp2.tsv
paste $OutDir/temp1.tsv $OutDir/temp2.tsv > $OutDir/experiment_all_gene_GO_annots_geneid.tsv
rm $OutDir/temp1.tsv
rm $OutDir/temp2.tsv

ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
AnnotTable=analysis/annotation_tables/WT_minion_gene_table.tsv
DEGs=analysis/Coexpression/Clust_Results_m70xq10z2qz49uwe/C0.txt
cat $DEGs | sed 's/.t.*//g' > $OutDir/C0_gene.txt
DEG2=analysis/Coexpression/enrichment/F.venenatum/WT_minion/whole_genome/C0/C0_gene.txt
AllGenes=$OutDir/Allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/score1.txt
Set2Genes=$OutDir/score2.txt
AllGenes=$OutDir/Cluster0.txt
cat $DEG2 | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes
$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt


OutDir=analysis/Coexpression/enrichment/F.venenatum/WT_minion/whole_genome/C1
mkdir -p $OutDir
InterProTSV=gene_pred/interproscan/F.venenatum/WT_minion/WT_minion_interproscan.tsv
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
$ProgDir/GO_table.py --interpro $InterProTSV > $OutDir/experiment_all_gene_GO_annots.tsv
cat $OutDir/experiment_all_gene_GO_annots.tsv | sed 's/.t.*//g' > $OutDir/temp1.tsv
cat $OutDir/experiment_all_gene_GO_annots.tsv | cut -f2 > $OutDir/temp2.tsv
paste $OutDir/temp1.tsv $OutDir/temp2.tsv > $OutDir/experiment_all_gene_GO_annots_geneid.tsv
rm $OutDir/temp1.tsv
rm $OutDir/temp2.tsv

ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
AnnotTable=analysis/annotation_tables/WT_minion_gene_table.tsv
DEGs=analysis/Coexpression/Clust_Results_m70xq10z2qz49uwe/C1.txt
cat $DEGs | sed 's/.t.*//g' > $OutDir/C1_gene.txt
DEG2=analysis/Coexpression/enrichment/F.venenatum/WT_minion/whole_genome/C1/C1_gene.txt
AllGenes=$OutDir/Allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/score1.txt
Set2Genes=$OutDir/score2.txt
AllGenes=$OutDir/Cluster0.txt
cat $DEG2 | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes
$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

OutDir=analysis/Coexpression/enrichment/F.venenatum/WT_minion/whole_genome/C2
mkdir -p $OutDir
InterProTSV=gene_pred/interproscan/F.venenatum/WT_minion/WT_minion_interproscan.tsv
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
$ProgDir/GO_table.py --interpro $InterProTSV > $OutDir/experiment_all_gene_GO_annots.tsv
cat $OutDir/experiment_all_gene_GO_annots.tsv | sed 's/.t.*//g' > $OutDir/temp1.tsv
cat $OutDir/experiment_all_gene_GO_annots.tsv | cut -f2 > $OutDir/temp2.tsv
paste $OutDir/temp1.tsv $OutDir/temp2.tsv > $OutDir/experiment_all_gene_GO_annots_geneid.tsv
rm $OutDir/temp1.tsv
rm $OutDir/temp2.tsv

ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
AnnotTable=analysis/annotation_tables/WT_minion_gene_table.tsv
DEGs=analysis/Coexpression/Clust_Results_m70xq10z2qz49uwe/C2.txt
cat $DEGs | sed 's/.t.*//g' > $OutDir/C2_gene.txt
DEG2=analysis/Coexpression/enrichment/F.venenatum/WT_minion/whole_genome/C2/C2_gene.txt
AllGenes=$OutDir/Allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/score1.txt
Set2Genes=$OutDir/score2.txt
AllGenes=$OutDir/Cluster0.txt
cat $DEG2 | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes
$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

OutDir=analysis/Coexpression/enrichment/F.venenatum/WT_minion/whole_genome/C3
mkdir -p $OutDir
InterProTSV=gene_pred/interproscan/F.venenatum/WT_minion/WT_minion_interproscan.tsv
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
$ProgDir/GO_table.py --interpro $InterProTSV > $OutDir/experiment_all_gene_GO_annots.tsv
cat $OutDir/experiment_all_gene_GO_annots.tsv | sed 's/.t.*//g' > $OutDir/temp1.tsv
cat $OutDir/experiment_all_gene_GO_annots.tsv | cut -f2 > $OutDir/temp2.tsv
paste $OutDir/temp1.tsv $OutDir/temp2.tsv > $OutDir/experiment_all_gene_GO_annots_geneid.tsv
rm $OutDir/temp1.tsv
rm $OutDir/temp2.tsv

ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
AnnotTable=analysis/annotation_tables/WT_minion_gene_table.tsv
DEGs=analysis/Coexpression/Clust_Results_m70xq10z2qz49uwe/C3.txt
cat $DEGs | sed 's/.t.*//g' > $OutDir/C3_gene.txt
DEG2=analysis/Coexpression/enrichment/F.venenatum/WT_minion/whole_genome/C3/C3_gene.txt
AllGenes=$OutDir/Allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/score1.txt
Set2Genes=$OutDir/score2.txt
AllGenes=$OutDir/Cluster0.txt
cat $DEG2 | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes
$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt







ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
AnnotTable=analysis/annotation_tables/WT_minion_gene_table.tsv
DEGs=analysis/Coexpression/Clust_Results_m70xq10z2qz49uwe/C1.txt
OutDir2=analysis/Coexpression/enrichment/F.venenatum/WT_minion/whole_genome/C1
mkdir -p $OutDir2
AllGenes=$OutDir2/Allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir2/score1.txt
Set2Genes=$OutDir2/score2.txt
AllGenes=$OutDir2/Cluster0.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes
$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir2 > $OutDir2/output.txt

ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
AnnotTable=analysis/annotation_tables/WT_minion_gene_table.tsv
DEGs=analysis/Coexpression/Clust_Results_m70xq10z2qz49uwe/C2.txt
OutDir2=analysis/Coexpression/enrichment/F.venenatum/WT_minion/whole_genome/C2
mkdir -p $OutDir2
AllGenes=$OutDir2/Allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir2/score1.txt
Set2Genes=$OutDir2/score2.txt
AllGenes=$OutDir2/Cluster2.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes
$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir2 > $OutDir2/output.txt

ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
AnnotTable=analysis/annotation_tables/WT_minion_gene_table.tsv
DEGs=analysis/Coexpression/Clust_Results_m70xq10z2qz49uwe/C3.txt
OutDir2=analysis/Coexpression/enrichment/F.venenatum/WT_minion/whole_genome/C3
mkdir -p $OutDir2
AllGenes=$OutDir2/Allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir2/score1.txt
Set2Genes=$OutDir2/score2.txt
AllGenes=$OutDir2/Cluster3.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes
$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir2 > $OutDir2/output.txt

WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_ALL/Wc1/Wc1vsWT53_D/DOWN
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL/Wc1vsWT53_D_down_LFC_names.txt
AllGenes=$OutDir/Wc1vsWT53_D_down_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Wc1vsWT53_D_down_DEGs.txt
Set2Genes=$OutDir/Wc1vsWT53_D_down2.txt
AllGenes=$OutDir/Wc1vsWT53_D_down.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```


grep -w -F -f ../../../../Clust_Results_m70xq10z2qz49uwe/C0.txt ../../../../../../analysis/annotation_tables/WT_minion_gene_table.tsv > C0/C0_annot.txt

grep -w -F -f ../../../../Clust_Results_m70xq10z2qz49uwe/C1.txt ../../../../../../analysis/annotation_tables/WT_minion_gene_table.tsv > C1/C1_annot.txt

grep -w -F -f ../../../../Clust_Results_m70xq10z2qz49uwe/C2.txt ../../../../../../analysis/annotation_tables/WT_minion_gene_table.tsv > C2/C2_annot.txt

grep -w -F -f ../../../../Clust_Results_m70xq10z2qz49uwe/C3.txt ../../../../../../analysis/annotation_tables/WT_minion_gene_table.tsv > C3/C3_annot.txt