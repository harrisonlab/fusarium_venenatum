# Generating an TSV file with sequencing information

```bash
#Antismash output correction
cat /projects/fusarium_venenatum/analysis/secondary_metabolites/antismash/F.venenatum/WT_minion_vAG/WT_antismash_results_secmet_genes.tsv | sed 's/;//p' | sed 's/;.*//p' | sed 's/Kin.*//p' > WT_antismash_results_secmet_genes_corrected.tsv


for GeneGff in $(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.gff3); do
    Strain=WT_minion
    Organism=F.venenatum
    Assembly=$(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
    TFs=$(ls analysis/transcription_factors/F.venenatum/WT_minion/WT_minion_TF_domains.tsv)
    InterPro=$(ls gene_pred/interproscan/F.venenatum/WT_minion/WT_minion_interproscan.tsv)
    Antismash=$(ls analysis/secondary_metabolites/antismash/F.venenatum/WT_minion_vAG/WT_antismash_results_secmet_genes_corrected.tsv)
    SwissProt=$(ls gene_pred/swissprot/F.venenatum/WT_minion/swissprot_vJun2020_tophit_parsed.tbl)
    OutDir=analysis/annotation_tables/$Organism/$Strain
    mkdir -p $OutDir
    GeneFasta=$(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta)
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Annotation_tables
    $ProgDir/build_annot_v1.py --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --TFs $TFs --InterPro $InterPro --Antismash $Antismash --Swissprot $SwissProt > $OutDir/"$Strain"_noDEGs_gene_table.tsv
done


```

# Cassis on all genes 

```bash
# screen -a
# ssh compute02
# srun -p medium --mem=4gb --pty bash

conda activate meme_suite

ProjDir=$(ls -d /projects/fusarium_venenatum)

WorkDir=$ProjDir/analysis/promoters/cassis/all_genes
mkdir -p $WorkDir
cd $WorkDir

AnnotTab=$(ls $ProjDir/analysis/annotation_tables/F.venenatum/WT_minion/WT_minion_noDEGs_gene_table.tsv)
Assembly=$(ls $ProjDir/repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/medaka_contigs_softmasked_repeatmasker_TPSI_appended.fa)
Genes=$(ls $ProjDir/gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.gff3)
Interpro=$(ls $ProjDir/gene_pred/interproscan/F.venenatum/WT_minion/WT_minion_interproscan.tsv)

cat $Genes | grep 'mRNA' | sed 's/ID=//g' | sed "s/;.*//g" | awk '{ print $9 "\t" $1 "\t" $4 "\t" $5 "\t" $7}' > cassis.tsv

CassisTSV=cassis.tsv




for Cluster in $(cat $AnnotTab | cut -f13 | grep 'contig' | sort -n -k3 -t'_' | sed 's/;.*//p' | uniq); do
echo $Cluster
mkdir $WorkDir/$Cluster
cat $AnnotTab | cut -f1,13 | grep -w "$Cluster" | cut -f1 | grep '.t1' > $WorkDir/$Cluster/headers.txt
for GeneID in $(cat $WorkDir/*/headers.txt); do
echo $GeneID
mkdir -p $WorkDir/$Cluster/$GeneID
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Promoter_analysis
OutDir=$WorkDir/$Cluster
Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
while [ $Jobs -gt 60 ]; do
sleep 5m
printf "."
Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
done		
printf "\n"
sbatch $ProgDir/cassis.sh $Assembly $CassisTSV $GeneID $OutDir
done
 
done
done
```

  cassis \
     --annotation cassis.tsv \
     --genome $Assembly \
     --anchor g13451.t1 \
     --dir contig_4_Cluster_9/g13459.t1 \
     --mismatches 0 \
     -v \
     --prediction \
     --num-cpus 30 \
     | tee 2>&1 $WorkDir/$Cluster/$GeneID/${GeneID}_log.txt