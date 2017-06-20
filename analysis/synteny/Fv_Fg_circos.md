## Preparing genome and .conf files

```bash
OutDir=analysis/circos/F.venenatum/Fv_Fg
mkdir -p $OutDir

Fv_genome=repeat_masked/F.venenatum/WT_ncbi/ncbi_submission/polished_contigs_unmasked.fa
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $Fv_genome --contig_prefix "A3_5_" > $OutDir/Fv_genome.txt

Fg_genome=assembly/external_group/F.graminearum/PH1/dna/Fusarium_graminearum.RR.dna.toplevel.fa
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $Fg_genome --contig_prefix "PH1_" > $OutDir/Fg_genome.txt

  cat $OutDir/Fv_genome.txt > $OutDir/Fv_Fg_genome.txt
  tac $OutDir/Fg_genome.txt >> $OutDir/Fv_Fg_genome.txt

  # COntigs smaller than 10Kb were removed
  cat $OutDir/Fv_Fg_genome.txt | grep -v -e 'PH1_Mt' -e 'PH1_HG970330' | grep -v -e "A3_5__contig_9" -e "A3_5_contig_10" -e "A3_5_contig_11" > $OutDir/Fv_Fg_genome_edited.txt
```
<!--
The order of contigs was changed manually using nano
```bash
cp $OutDir/Fv_Fg_genome_edited.txt $OutDir/Fv_Fg_genome_edited2.txt
nano $OutDir/Fv_Fg_genome_edited2.txt
cp $OutDir/Fv_Fg_genome_edited2.txt $OutDir/Fv_Fg_genome_final.txt
``` -->


```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/orthology2circos_ribbons.py --orthology analysis/orthology/orthomcl/Fv_vs_Fg/Fv_vs_Fg_orthogroups.txt --name1 A3_5 --gff1 gene_pred/final/F.venenatum/WT_ncbi/final/final_genes_appended.gff3 --name2 PH1 --gff2 assembly/external_group/F.graminearum/PH1/gff/Fusarium_graminearum.RR.34.gff3 \
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
  cat $OutDir/Fv_genome.txt | grep -v -e 'Fg_Mt' -e 'Fg_HG970330' | grep -v -e "Fv_contig_9" -e "Fv_contig_10" -e "Fv_contig_11" | grep -f $OutDir/Fv_syntenous_contigs.txt | cut -f6 -d ' ' | awk '{s+=$1} END {print s}'
```

Contig order was selected by taking the first line of that file and then also
taking the reversed order of Fg contigs using the command:

```bash
cat $OutDir/Fv_Fg_contig_orientation.txt | grep -A1 'Order of all seen contigs'
cat $OutDir/Fv_Fg_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" | tr -d '\n' | sed 's/A3_5/, A3_5/g'  > tmp.txt
cat $OutDir/Fv_Fg_genome_edited.txt | grep 'A3_5' | grep -v -e 'Fg_Mt' -e 'Fg_HG970330' | grep -v -e "Fv_contig_9" -e "Fv_contig_10" -e "Fv_contig_11" | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/Fv/, Fv/g'
cat $OutDir/Fv_Fg_genome_edited.txt | grep 'PH1' | cut -f3 -d ' ' | tr -d '\n' | sed 's/PH1/, PH1/g' >> tmp.txt
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
