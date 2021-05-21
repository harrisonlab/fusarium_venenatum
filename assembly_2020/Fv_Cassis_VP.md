# Generating an TSV file with sequencing information

```bash
#Antismash output correction
cat analysis/secondary_metabolites/antismash/F.venenatum/WT_minion_VP/WT_antismash_results_secmet_genes.tsv | sed 's/;//p' | sed 's/;.*//p' | sed 's/Kin.*//p' > analysis/secondary_metabolites/antismash/F.venenatum/WT_minion_VP/WT_antismash_results_secmet_genes_corrected.tsv


for GeneGff in $(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.gff3); do
Strain=WT_minion
Organism=F.venenatum
Assembly=$(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
TFs=$(ls analysis/transcription_factors/F.venenatum/WT_minion/WT_minion_TF_domains.tsv)
InterPro=$(ls gene_pred/interproscan/F.venenatum/WT_minion/WT_minion_interproscan.tsv)
Antismash=$(ls analysis/secondary_metabolites/antismash/F.venenatum/WT_minion_VP/WT_antismash_results_secmet_genes_corrected.tsv)
#Smurf=$(ls analysis/secondary_metabolites/smurf/F.venenatum/WT_minion/WT_minion_smurf_secmet_genes.tsv) # I added cassis genes manually
SwissProt=$(ls gene_pred/swissprot/F.venenatum/WT_minion/swissprot_vJun2020_tophit_parsed.tbl)
OutDir=analysis/annotation_tables_VP/$Organism/$Strain
mkdir -p $OutDir
GeneFasta=$(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta)
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Annotation_tables
$ProgDir/build_annot_v1.py --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --TFs $TFs --InterPro $InterPro --Antismash $Antismash --Swissprot $SwissProt > $OutDir/"$Strain"_noDEGs_gene_table.tsv
done
```
