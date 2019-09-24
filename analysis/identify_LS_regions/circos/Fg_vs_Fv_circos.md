```bash
  OutDir=analysis/genome_alignment/circos/Fg_vs_Fv_cov_plot
  mkdir -p $OutDir

  # Ref_genome=$(ls ../../../../../home/groups/harrisonlab/project_files/fusarium_venenatum/repeat_masked/F.venenatum/WT_minion/minion_submission/*_contigs_unmasked.fa)
  Ref_genome=$(ls /home/groups/harrisonlab/project_files/fusarium_venenatum/repeat_masked/F.venenatum/WT/illumina_assembly_ncbi/*_contigs_unmasked.fa)
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/fasta2circos.py --genome $Ref_genome --contig_prefix "" | grep -P "\d\d\d\d\d chr" > $OutDir/ref_genome.txt
```


Calculate coverage over 10kb windows.
Parse read depth files:

```bash
OutDir=analysis/genome_alignment/circos/Fg_vs_Fv_cov_plot
for File in $(ls analysis/genome_alignment/bwa/*/*/vs_WT/*_vs_WT_depth_10kb.tsv); do
  Strain=$(echo $File | cut -f5 -d '/')
  echo $Strain
  cat $File | awk '{print $1,$2-1000,$2,$3,$4}' OFS='\t' | cut -f1,2,3,4 > $OutDir/${Strain}_vs_ref_unmasked_scatterplot.tsv
done
```

```bash
OutDir=analysis/genome_alignment/circos/Fg_vs_Fv_cov_plot
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_venenatum/analysis/identify_LS_regions/circos
circos -conf $ProgDir/Fg_vs_Fv_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/Fg_vs_Fv_cov_unmasked_circos.png
mv $OutDir/circos.svg $OutDir/Fg_vs_Fv_cov_unmasked_circos.svg
# mv $OutDir/circos.png $OutDir/Y-11545_v2_cov_masked_circos.png
# mv $OutDir/circos.svg $OutDir/Y-11545_v2_cov_masked_circos.svg
ls $PWD/$OutDir/Fg_vs_Fv*_circos.png
```
