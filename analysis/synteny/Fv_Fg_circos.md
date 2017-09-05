## Preparing genome and .conf files

```bash
OutDir=analysis/circos/F.venenatum/Fv_Fg
mkdir -p $OutDir

Fv_genome=$(ls repeat_masked/F.venenatum/WT/illumina_assembly_ncbi/WT_contigs_unmasked.fa)
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $Fv_genome --contig_prefix "A3_5_" > $OutDir/Fv_genome.txt

Fg_genome=$(ls assembly/external_group/F.graminearum/PH1/dna/Fusarium_graminearum.RR1.dna.toplevel.fa)
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $Fg_genome --contig_prefix "PH1_" > $OutDir/Fg_genome.txt

  cat $OutDir/Fv_genome.txt > $OutDir/Fv_Fg_genome.txt
  tac $OutDir/Fg_genome.txt >> $OutDir/Fv_Fg_genome.txt

  # Contigs smaller than 10Kb were removed
  cat $OutDir/Fv_Fg_genome.txt | grep -v -e 'PH1_Mt' -e 'PH1_HG970330' \
  | grep -v -e "A3_5_contig_87" \
  | grep -v -e "A3_5_contig_88" \
  | grep -v -e "A3_5_contig_89" \
  | grep -v -e "A3_5_contig_90" \
  | grep -v -e "A3_5_contig_91" \
  | grep -v -e "A3_5_contig_92" \
  | grep -v -e "A3_5_contig_93" \
  | grep -v -e "A3_5_contig_94" \
  | grep -v -e "A3_5_contig_95" \
  | grep -v -e "A3_5_contig_96" \
  | grep -v -e "A3_5_contig_97" \
  | grep -v -e "A3_5_contig_98" \
  | grep -v -e "A3_5_contig_99" \
  | grep -v -e "A3_5_contig_100" \
  | grep -v -e "A3_5_contig_101" \
  | grep -v -e "A3_5_contig_102" \
  | grep -v -e "A3_5_contig_103" \
  | grep -v -e "A3_5_contig_104" \
  | grep -v -e "A3_5_contig_105" \
  > $OutDir/Fv_Fg_genome_edited.txt
```
<!--
The order of contigs was changed manually using nano
```bash
cp $OutDir/Fv_Fg_genome_edited.txt $OutDir/Fv_Fg_genome_edited2.txt
nano $OutDir/Fv_Fg_genome_edited2.txt
cp $OutDir/Fv_Fg_genome_edited2.txt $OutDir/Fv_Fg_genome_final.txt
``` -->


```bash
  Orthology=$(ls analysis/orthology/orthomcl/Fv_vs_Fg/Fv_vs_Fg_orthogroups.txt)
  Gff1=$(ls gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.gff3)
  Gff2=$(ls assembly/external_group/F.graminearum/PH1/gff/Fusarium_graminearum.RR1.36.gff3)
  cat $Gff2 | grep -v '#' > $OutDir/tmp.gff3
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/orthology2circos_ribbons.py --orthology $Orthology --name1 A3_5 --gff1 $Gff1 --name2 PH1 --gff2 $OutDir/tmp.gff3 \
   | sort -k4,5 -V \
   > $OutDir/Fv_Fg_links.txt
  # Links to Fg LS contigs 3, 6, 14 and 15 were coloured black
  # cat $OutDir/Fv_Fg_links.txt \
  #   | sed '/4287_CM000591.1/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000594.1/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000602.2/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000603.1/ s/$/\tcolor=black/' \
  #   > $OutDir/Fv_Fg_links_edited.txt
  cat $OutDir/Fv_Fg_links.txt > $OutDir/Fv_Fg_links_edited.txt
```

A file showing contig orientations was made:
```bash
cat $OutDir/Fv_Fg_links_edited.txt | cut -f1 | uniq > $OutDir/Fv_contig_order.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/find_contig_orientation.py --links_file $OutDir/Fv_Fg_links_edited.txt > $OutDir/Fv_Fg_contig_orientation.txt
```

The number of bp in syntenous contigs was identified using:

```bash
  cat $OutDir/Fv_Fg_contig_orientation.txt | tail -n3 | grep -v 'orientation' | sed 's/, /\n/g' > $OutDir/Fv_syntenous_contigs.txt
  cat $OutDir/Fv_genome.txt | grep -v -e 'PH1_Mt' -e 'PH1_HG970330' \
  | grep -v -e "A3_5_contig_87" \
  | grep -v -e "A3_5_contig_88" \
  | grep -v -e "A3_5_contig_89" \
  | grep -v -e "A3_5_contig_90" \
  | grep -v -e "A3_5_contig_91" \
  | grep -v -e "A3_5_contig_92" \
  | grep -v -e "A3_5_contig_93" \
  | grep -v -e "A3_5_contig_94" \
  | grep -v -e "A3_5_contig_95" \
  | grep -v -e "A3_5_contig_96" \
  | grep -v -e "A3_5_contig_97" \
  | grep -v -e "A3_5_contig_98" \
  | grep -v -e "A3_5_contig_99" \
  | grep -v -e "A3_5_contig_100" \
  | grep -v -e "A3_5_contig_101" \
  | grep -v -e "A3_5_contig_102" \
  | grep -v -e "A3_5_contig_103" \
  | grep -v -e "A3_5_contig_104" \
  | grep -v -e "A3_5_contig_105" \
  | grep -f $OutDir/Fv_syntenous_contigs.txt | cut -f6 -d ' ' | awk '{s+=$1} END {print s}'
```

```
37597212
```

Contig order was selected by taking the first line of that file and then also
taking the reversed order of Fg contigs using the command:

```bash
cat $OutDir/Fv_Fg_contig_orientation.txt | grep -A1 'Order of all seen contigs'
cat $OutDir/Fv_Fg_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" | tr -d '\n' | sed 's/A3_5/, A3_5/g'  > tmp.txt
# cat $OutDir/Fv_Fg_genome_edited.txt | grep 'A3_5' | grep -v -e 'PH1_Mt' -e 'PH1_HG970330' \
# | grep -v -e "A3_5_contig_87" \
# | grep -v -e "A3_5_contig_88" \
# | grep -v -e "A3_5_contig_89" \
# | grep -v -e "A3_5_contig_90" \
# | grep -v -e "A3_5_contig_91" \
# | grep -v -e "A3_5_contig_92" \
# | grep -v -e "A3_5_contig_93" \
# | grep -v -e "A3_5_contig_94" \
# | grep -v -e "A3_5_contig_95" \
# | grep -v -e "A3_5_contig_96" \
# | grep -v -e "A3_5_contig_97" \
# | grep -v -e "A3_5_contig_98" \
# | grep -v -e "A3_5_contig_99" \
# | grep -v -e "A3_5_contig_100" \
# | grep -v -e "A3_5_contig_101" \
# | grep -v -e "A3_5_contig_102" \
# | grep -v -e "A3_5_contig_103" \
# | grep -v -e "A3_5_contig_104" \
# | grep -v -e "A3_5_contig_105" \
# | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/A3_5/, A3_5/g'
# echo ""
cat $OutDir/Fv_Fg_genome_edited.txt | grep 'PH1' | cut -f3 -d ' ' | tr -d '\n' | sed 's/PH1/, PH1/g' >> tmp.txt
```

```
A3_5_contig_18, A3_5_contig_30, A3_5_contig_44, A3_5_contig_29, A3_5_contig_46, A3_5_contig_43, A3_5_contig_70, A3_5_contig_60, A3_5_contig_54, A3_5_contig_7, A3_5_contig_76, A3_5_contig_58, A3_5_contig_52, A3_5_contig_68, A3_5_contig_84, A3_5_contig_16, A3_5_contig_28, A3_5_contig_67, A3_5_contig_21, A3_5_contig_10, A3_5_contig_25, A3_5_contig_37, A3_5_contig_65, A3_5_contig_77, A3_5_contig_74, A3_5_contig_82, A3_5_contig_26, A3_5_contig_75, A3_5_contig_78, A3_5_contig_34, A3_5_contig_45, A3_5_contig_53, A3_5_contig_14, A3_5_contig_32, A3_5_contig_24, A3_5_contig_17, A3_5_contig_35, A3_5_contig_5, A3_5_contig_50, A3_5_contig_69, A3_5_contig_15, A3_5_contig_31, A3_5_contig_6, A3_5_contig_86, A3_5_contig_38, A3_5_contig_81, A3_5_contig_61, A3_5_contig_62, A3_5_contig_12, A3_5_contig_1, A3_5_contig_33, A3_5_contig_22, A3_5_contig_66, A3_5_contig_57, A3_5_contig_13, A3_5_contig_9, A3_5_contig_59, A3_5_contig_47, A3_5_contig_23, A3_5_contig_19, A3_5_contig_64, A3_5_contig_79, A3_5_contig_51, A3_5_contig_20, A3_5_contig_41, A3_5_contig_49, A3_5_contig_63, A3_5_contig_2, A3_5_contig_55, A3_5_contig_80, A3_5_contig_11, A3_5_contig_48, A3_5_contig_36, A3_5_contig_42, A3_5_contig_40, A3_5_contig_83, A3_5_contig_39, A3_5_contig_27, A3_5_contig_4, A3_5_contig_73, A3_5_contig_3, A3_5_contig_8, A3_5_contig_56, A3_5_contig_71, A3_5_contig_72, A3_5_contig_85, PH1_4, PH1_3, PH1_2, PH1_1

```

Contig orientation was used to edit the circos .conf file manually
<!--
## Preparing Effector plots

# Plot location of Fus2 mimps and secreted effectorP genes as a scatterplot

```bash
  GffMimp=analysis/mimps/F.venenatum/WT_ncbi/WT_ncbi_mimps.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffMimp --feature MIMP_motif --value 1 | sed -e 's/^/Fv_/g' > $OutDir/Fv_mimp_plot.txt
  GffCAZY=gene_pred/CAZY/F.venenatum/WT_ncbi/WT_ncbi_CAZY_secreted.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffCAZY --feature gene --value 1 | sed -e 's/^/Fv_/g' > $OutDir/Fv_CAZY_plot.txt
  GffEfFv=analysis/effectorP/F.venenatum/WT_ncbi/F.venenatum_WT_ncbi_EffectorP_secreted.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffEfFv --feature gene --value 1 | sed -e 's/^/Fv_/g' > $OutDir/Fv_effectorP_plot.txt
  GffAntiSmash=analysis/antismash/F.venenatum/WT_ncbi/WT_ncbi_secondary_metabolite_regions.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffAntiSmash --feature indole indole-nrps nrps nrps-t1pks other t1pks t1pks-nrps t3pks terpene --value 1 | sed -e 's/^/Fv_/g' > $OutDir/Fv_antismash_plot.txt

  BlastHits=analysis/blast_homology/F.venenatum/WT_ncbi/WT_ncbi_Fo_path_genes_CRX.fa_homologs.gff
  GffSix=$OutDir/Fv_SIX.gff
  cat $BlastHits | grep -v -e 'MIMP' -e 'C5' -e 'CRX' > $GffSix
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffSix --feature SIX_homolog --value 1 | sed -e 's/^/Fv_/g' > $OutDir/Fv_SIX_plot.txt
``` -->

## Running circos



```bash
Conf=$(ls /home/armita/git_repos/emr_repos/scripts/fusarium_venenatum/analysis/synteny/Fv_Fg_circos.conf)
circos -conf $Conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/Fv_Fg_circos.png
mv $OutDir/circos.svg $OutDir/Fv_Fg_circos.svg
```


# Further analysis of non-syntenous regions regions

The number of MIMPs and effectors in LS regions were identified:


```bash
cat $OutDir/Fv_mimp_plot.txt | grep -v -f $OutDir/Fv_syntenous_contigs.txt | wc -l

```
