# Fusarium venenatum Genome Assembly

======================================

This document details the commands used to assemble and annotate a reference genome of fusarium venenatum

The following includes new approaches in the assembly of Fvenenatum genomes.

Note - all this work was performed in the directory: /projects/fusarium_venenatum

## Data management

Data from two sequencing runs were merged together

```bash
mkdir -p qc_dna/minion/F.venenatum/WT

cat F.venenatum_WT_07-03-17_albacore_v2.02_trim.fastq.gz F.venenatum_WT_18-07-17_albacore_v2.02_trim.fastq.gz > /projects/fusarium_venenatum/qc_dna/minion/F.venenatum/WT/WT_minion_allfiles.fastq.gz
```


## Miniasm

```bash
# conda packages required: miniasm, minimap2, bbmap, busco
conda activate olc_assemblers
```

## Assembly of uncorrected reads. Racon required after.

```bash
cd qc_dna/minion/F.venenatum/WT
# If reads have same names or same part splitted by space, fix them using rename.sh from bbtools.
rename.sh in=WT_minion_allfiles.fastq.gz out=WT_minion_allfiles.fasta prefix=WT_minion

# Fast all-against-all overlap of raw reads
# Overlap for MinION reads (or use "-x ava-pb" for Pacbio read overlapping)
minimap2 -x ava-ont -t8 WT_minion_allfiles.fasta WT_minion_allfiles.fasta | gzip -1 > WT_minion_fastq_allfiles.paf.gz

# Concatenate pieces of read sequences to generate the final sequences.
# Thus the per-base error rate is similar to the raw input reads. Make sure you correct your reads.
# Layout
miniasm -f WT_minion_allfiles.fasta WT_minion_fastq_allfiles.paf.gz > reads.gfa

# Convert gfa file to fasta file.
awk '/^S/{print ">"$2"\n"$3}' reads.gfa | fold > assembly/miniasm/F.venenatum/WT_minion/WT_minion.fa
```

Quast, busco and kat were run to assess the assembly quality (optionally at this stage, racon will correct errors after)

```bash
# Python 2.7 is needed to install Quast
  ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/assembly_qc
  for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/WT_minion.fa); do
    OutDir=$(dirname $Assembly)
    sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
mv assembly/miniasm/F.venenatum/WT_minion/WT_minion.fa assembly/miniasm/F.venenatum/WT_minion/WT_minion_miniasm.fa
  for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/WT_minion_miniasm.fa); do
    Strain=WT_minion
    Organism=F.venenatum
    echo "$Organism - $Strain"
    ProgDir=/home/gomeza/git_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /projects/oldhome/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd9
    sbatch $ProgDir/sub_busco.sh $Assembly $BuscoDB $OutDir
  done
```




Kat tool was used to assess level of error and duplications in the genome assemblies generated. Mapleson et al., 2016.

```bash
#Editing....
for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/WT_minion.fa); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev )
  Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
  echo "$Organism - $Strain"
  IlluminaDir=$(ls -d qc_dna/paired/*/Hg199)
  cat $IlluminaDir/F/*_trim.fq.gz > $IlluminaDir/F/F_trim_appended.fq.gz
  cat $IlluminaDir/R/*_trim.fq.gz > $IlluminaDir/R/R_trim_appended.fq.gz
  ReadsF=$(ls $IlluminaDir/F/F_trim_appended.fq.gz)
  ReadsR=$(ls $IlluminaDir/R/R_trim_appended.fq.gz)
  OutDir=$(dirname $Assembly)/kat
  Prefix="${Strain}"
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/kat
  qsub $ProgDir/sub_kat.sh $Assembly $ReadsF $ReadsR $OutDir $Prefix 200
done
```

after KAT jobs have finished running, then remove appended trimmed reads
```bash
rm qc_dna/paired/*/*/*/F_trim_appended.fq.gz
rm qc_dna/paired/*/*/*/R_trim_appended.fq.gz
```



## Error correction using racon:

```bash
  for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/WT_minion.fa); do
    Strain=WT
    Organism=F.venenatum
    echo "$Organism - $Strain"
    ReadsFq=$(ls qc_dna/minion/F.venenatum/WT/WT_minion_allfiles.fastq.gz)
    Iterations=10
    OutDir=$(dirname $Assembly)"/racon_$Iterations"
    ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/racon
    sbatch $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
  done
```

Quast and busco were run to assess the effects of racon on assembly quality:


```bash
# Python 2.7 is needed to install Quast
  ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/assembly_qc
  for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/racon_10/WT_minion_racon_round_10.fasta); do
    OutDir=$(dirname $Assembly)
    sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```


```bash
ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
touch tmp.txt
for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/racon_10/WT_minion_racon_round_10.fasta); do
  OutDir=$(dirname $Assembly)
  $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/WT_miniasm_racon10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
rm tmp.txt
```

```bash
  for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/racon_10/WT_miniasm_racon10_renamed.fasta); do
    Strain=WT_minion
    Organism=F.venenatum
    echo "$Organism - $Strain"
    ProgDir=/home/gomeza/git_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /projects/oldhome/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd9
    sbatch $ProgDir/sub_busco.sh $Assembly $BuscoDB $OutDir
  done
```

Assemblies were polished using Pilon
Although three libraries were available, the first contained a relatively small amount of data and was not used for correction.


```bash
for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/racon_10/WT_miniasm_racon10_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev | sed 's/_minion//g')
IlluminaDir=$(ls -d ../oldhome/groups/harrisonlab/project_files/fusarium_venenatum/qc_dna/paired/$Organism/$Strain)
echo $Strain
echo $Organism
TrimF2_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n2 | tail -n1);
TrimR2_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n2 | tail -n1);
TrimF3_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n3 | tail -n1);
TrimR3_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n3 | tail -n1);
echo $TrimF2_Read
echo $TrimR2_Read
echo $TrimF3_Read
echo $TrimR3_Read
OutDir=$(dirname $Assembly)
Iterations=10
ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/pilon
sbatch $ProgDir/sub_pilon_2_libs.sh $Assembly $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir $Iterations
done
```

```bash
# Python 2.7 is needed to install Quast
  ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/assembly_qc
  for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/racon_10/pilon_5.fasta); do
    OutDir=$(dirname $Assembly)/pilon5
    sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```


```bash
  for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/racon_10/pilon_5.fasta); do
    Strain=WT_minion
    Organism=F.venenatum
    echo "$Organism - $Strain"
    ProgDir=/home/gomeza/git_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /projects/oldhome/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd9
    sbatch $ProgDir/sub_busco.sh $Assembly $BuscoDB $OutDir
  done
```



## SMARTdenovo

```bash
  for TrimReads in $(ls qc_dna/minion/F.venenatum/WT/WT_minion_allfiles.fastq.gz); do
    Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
    Strain=WT
    Prefix="$Strain"_smartdenovo
    OutDir=assembly/SMARTdenovo/F.venenatum/"$Strain"
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/SMARTdenovo
    sbatch $ProgDir/sub_SMARTdenovo.sh $TrimReads $Prefix $OutDir
  done
```

Quast, busco and kat were run to assess the assembly quality.

```bash
# Python 2.7 is needed to install Quast
  ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/assembly_qc
  for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT/WT_smartdenovo.dmo.lay.utg); do
    OutDir=$(dirname $Assembly)
    sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT/WT_smartdenovo.dmo.lay.utg); do
    ReadsFq=$(ls qc_dna/minion/F.venenatum/WT/WT_minion_allfiles.fastq.gz)
    Iterations=10
    OutDir=$(dirname $Assembly)"/racon_$Iterations"
    ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/racon
    sbatch $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
  done
```

```bash
ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
touch tmp.txt
for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_minion/racon_10/WT_smartdenovo_racon_round_10.fasta); do
  OutDir=$(dirname $Assembly)
  $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/WT_minion_racon10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
rm tmp.txt
```

```bash
# Python 2.7 is needed to install Quast
  ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/assembly_qc
  for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_minion/racon_10/WT_minion_racon10_renamed.fasta); do
    OutDir=$(dirname $Assembly)
    sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_minion/racon_10/WT_minion_racon10_renamed.fasta); do
    Strain=WT_minion
    Organism=F.venenatum
    echo "$Organism - $Strain"
    ProgDir=/home/gomeza/git_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /projects/oldhome/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd9
    sbatch $ProgDir/sub_busco.sh $Assembly $BuscoDB $OutDir
  done
```



Assemblies were polished using Pilon
Although three libraries were available, the first contained a relatively small amount of data and was not used for correction.


```bash
for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_minion/racon_10/WT_minion_racon10_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev | sed 's/_minion//g')
IlluminaDir=$(ls -d ../oldhome/groups/harrisonlab/project_files/fusarium_venenatum/qc_dna/paired/$Organism/$Strain)
echo $Strain
echo $Organism
TrimF2_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n2 | tail -n1);
TrimR2_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n2 | tail -n1);
TrimF3_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n3 | tail -n1);
TrimR3_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n3 | tail -n1);
echo $TrimF2_Read
echo $TrimR2_Read
echo $TrimF3_Read
echo $TrimR3_Read
OutDir=$(dirname $Assembly)
Iterations=10
ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/pilon
sbatch $ProgDir/sub_pilon_2_libs.sh $Assembly $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir $Iterations
done
```

```bash
ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
touch tmp.txt
for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/racon_10/pilon_10.fasta); do
  OutDir=$(dirname $Assembly)
  $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/WT_miniasm_pilon10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
rm tmp.txt
```


```bash
# Python 2.7 is needed to install Quast
  ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/assembly_qc
  for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/racon_10/pilon/WT_miniasm_pilon10_renamed.fasta); do
    OutDir=$(dirname $Assembly)/pilon6
    sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```


```bash
  for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/racon_10/pilon/WT_miniasm_pilon10_renamed.fasta); do
    #Strain=WT_minion
    #Organism=F.venenatum
    #echo "$Organism - $Strain"
    ProgDir=/home/gomeza/git_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /projects/oldhome/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd9
    sbatch $ProgDir/sub_busco.sh $Assembly $BuscoDB $OutDir
  done
```
```bash
ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
touch tmp.txt
for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_minion/racon_10/pilon/pilon_10.fasta); do
  OutDir=$(dirname $Assembly)
  $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/WT_SMART_pilon10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
rm tmp.txt
```



```bash
# Python 2.7 is needed to install Quast
  ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/assembly_qc
  for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_minion/racon_10/pilon/WT_SMART_pilon10_renamed.fasta); do
    OutDir=$(dirname $Assembly)/pilon6
    sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```


```bash
  for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_minion/racon_10/pilon/WT_SMART_pilon10_renamed.fasta); do
    ProgDir=/home/gomeza/git_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /projects/oldhome/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd9
    sbatch $ProgDir/sub_busco.sh $Assembly $BuscoDB $OutDir
  done
```





```bash
for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_minion/racon_10/WT_minion_racon10_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ReadDir=../../home/gomeza/rawdata4nanopolish/F.venenatum/WT
#ScratchDir=/data/scratch/nanopore_tmp_data/Fven
#Fast5Dir1=$ScratchDir/F.venenatum_WT_07-03-17/workspace/pass
#Fast5Dir2=$ScratchDir/F.venenatum_WT_18-07-17/workspace/pass
#nanopolish index -d $Fast5Dir1 -d $Fast5Dir2 $ReadDir/"$Strain"_concatenated_reads_filtered.fastq
OutDir=nanopolish_smart
mkdir -p $OutDir
ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/nanopolish
sbatch $ProgDir/sub_bwa_nanopolish_himem.sh $Assembly $ReadDir/"$Strain"_concatenated_reads_filtered.fastq $OutDir/nanopolish
done
```



```bash

  for Assembly in $(ls assembly_vAG/canu_2step/canu_minion/N.ditissima/Hg199/pilon/*10.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```


## Assembly correction with nanopolish

```bash
ReadDir=rawdata4nanopolish/F.venenatum/WT
mkdir -p $ReadDir
Strain=WT_minion
ReadsFq1=$(ls /projects/oldhome/groups/harrisonlab/project_files/fusarium_venenatum/raw_dna/minion/F.venenatum/WT/F.venenatum_WT_07-03-17_albacore_v2.02.fastq.gz)
ReadsFq2=$(ls /projects/oldhome/groups/harrisonlab/project_files/fusarium_venenatum/raw_dna/minion/F.venenatum/WT/F.venenatum_WT_18-07-17_albacore_v2.02.fastq.gz)
cat $ReadsFq1 $ReadsFq2 | gunzip -cf > $ReadDir/"$Strain"_concatenated_reads.fastq

/home/gomeza/git_repos/tools/seq_tools/assemblers/nanopolish/nanopolish_remove_dup_reads.py --fastq $ReadDir/"$Strain"_concatenated_reads.fastq --out $ReadDir/"$Strain"_concatenated_reads_filtered.fastq
```
```bash
for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/racon_10/WT_minion_racon_round_10.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ReadDir=../../home/gomeza/rawdata4nanopolish/F.venenatum/WT
#ScratchDir=/data/scratch/nanopore_tmp_data/Fven
#Fast5Dir1=$ScratchDir/F.venenatum_WT_07-03-17/workspace/pass
#Fast5Dir2=$ScratchDir/F.venenatum_WT_18-07-17/workspace/pass
#nanopolish index -d $Fast5Dir1 -d $Fast5Dir2 $ReadDir/"$Strain"_concatenated_reads_filtered.fastq
OutDir=nanopolish_bwa
mkdir -p $OutDir
ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/nanopolish
sbatch $ProgDir/sub_bwa_nanopolish.sh $Assembly $ReadDir/"$Strain"_concatenated_reads_filtered.fastq $OutDir/nanopolish
done
```


Split the assembly into 50Kb fragments an submit each to the cluster for nanopolish correction

```bash
for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/racon_10/WT_miniasm_racon10_renamed.fasta ); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly/nanopolish_variants)
RawReads=$(ls ../../home/gomeza/rawdata4nanopolish/F.venenatum/WT_minion/"$Strain"_concatenated_reads_filtered.fastq)
AlignedReads=$(ls $OutDir/nanopolish/reads.sorted.bam)

#NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
#nanopolish_makerange.py $Assembly > $OutDir/nanopolish/nanopolish_range.txt

Ploidy=1
#echo "nanopolish log:" > nanopolish_log.txt
#for Region in $(cat $OutDir/nanopolish/nanopolish_range.txt | head -n1); do
#Jobs=$(squeue | grep 'nanopo' | grep 'qw' | wc -l)
#while [ $Jobs -gt 1 ]; do
#sleep 1m
#printf "."
#Jobs=$(squeue | grep 'nanopo' | grep 'qw' | wc -l)
#done		
#printf "\n"
#echo $Region
#echo $Region >> nanopolish_log.txt
ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/nanopolish
sbatch --array=1-1000%50 $ProgDir/sub_nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $Region $OutDir/$Region
done
done
```

## Canu

```bash
	for Reads in $(ls qc_dna/minion/F.venenatum/WT/WT_minion_allfiles.fastq.gz); do
  	GenomeSz="37m"
  	Strain=WT
  	Organism=$(echo $Reads | rev | cut -f3 -d '/' | rev)
  	Prefix="$Strain"_canu
  	OutDir=assembly/canu/$Organism/"$Strain"
  	ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/canu
  	sbatch $ProgDir/submit_canu_minion.sh $Reads $GenomeSz $Prefix $OutDir
  done
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
# Python 2.7 is needed to install Quast
  ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/assembly_qc
  for Assembly in $(ls assembly/canu/F.venenatum/WT/WT_canu.contigs.fasta); do
    OutDir=$(dirname $Assembly)
    sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly/canu/F.venenatum/WT_minion/WT_canu.contigs.fasta); do
    Strain=WT
    Organism=F.venenatum
    echo "$Organism - $Strain"
    ProgDir=/home/gomeza/git_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /projects/oldhome/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd9
    sbatch $ProgDir/sub_busco.sh $Assembly $BuscoDB $OutDir
  done
```






## Previous assemblies qc

Quast, busco and kat were run to assess the assembly quality.

```bash
# Python 2.7 is needed to install Quast
  ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/assembly_qc
  for Assembly in $(ls assembly/previous_versions/F.venenatum/*/*.fa); do
    OutDir=$(dirname $Assembly)
    sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly/previous_versions/F.venenatum/WT_minion/WT_albacore_v2_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=WT_minion
    Organism=F.venenatum
    echo "$Organism - $Strain"
    ProgDir=/home/gomeza/git_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /projects/oldhome/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd9
    sbatch $ProgDir/sub_busco.sh $Assembly $BuscoDB $OutDir
  done
```

```bash
  for Assembly in $(ls assembly/previous_versions/F.venenatum/WT/WT_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=WT_minion
    Organism=F.venenatum
    echo "$Organism - $Strain"
    ProgDir=/home/gomeza/git_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /projects/oldhome/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd9
    sbatch $ProgDir/sub_busco.sh $Assembly $BuscoDB $OutDir
  done
```

flye
  for TrimReads in $(ls qc_dna/minion/F.venenatum/WT/WT_minion_allfiles.fastq.gz); do
    Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
    Strain=WT_minion
    Prefix="$Strain"_flye
    OutDir=assembly/flye/F.venenatum/"$Strain"
    mkdir -p $OutDir
    Size=
    ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/SMARTdenovo
    sbatch $ProgDir/flye.sh $TrimReads $Prefix $OutDir $Size
  done

  ```bash
# Python 2.7 is needed to install Quast
ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/assembly_qc
for Assembly in $(ls flye/assembly.fasta); do
OutDir=$(dirname $Assembly)
sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls flye/assembly.fasta); do
Strain=WT_minion
Organism=F.venenatum
echo "$Organism - $Strain"
ProgDir=/home/gomeza/git_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /projects/oldhome/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd9
sbatch $ProgDir/sub_busco.sh $Assembly $BuscoDB $OutDir
done
```


## Assembly correction with nanopolish

```bash
  ReadDir=rawdata4nanopolish/$Organism/$Strain
  mkdir -p $ReadDir
  ReadsFq1=$(ls /path/to/raw/basecalled/minion/reads/e.g./F.venenatum_WT_07-03-17_albacore_v2.02.fastq.gz)
  ReadsFq2=$(ls /path/to/raw/basecalled/minion/reads/e.g./F.venenatum_WT_18-07-17_albacore_v2.02.fastq.gz)
  cat $ReadsFq1 $ReadsFq2 | gunzip -cf > $ReadDir/"$Strain"_concatenated_reads.fastq

# Remove duplicate reads
  /home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assembler/nanopolish_remove_dup_reads.py --fastq $ReadDir/"$Strain"_concatenated_reads.fastq --out $ReadDir/"$Strain"_concatenated_reads_filtered.fastq

# Build an index mapping from basecalled reads to the signals measured by the sequencer
  ScratchDir=/data/scratch/nanopore_tmp_data/$Organism/$Strain
  Fast5Dir1=$ScratchDir/path/to/Fast5Files/workspace/pass
  Fast5Dir2=$ScratchDir/path/to/Fast5Files/workspace/pass
  nanopolish index -d $Fast5Dir1 -d $Fast5Dir2 $ReadDir/"$Strain"_concatenated_reads_filtered.fastq
```

```bash
  for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/racon_10/WT_minion_racon_round_10.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    ReadDir=../../home/gomeza/rawdata4nanopolish/F.venenatum/WT
    OutDir=nanopolish_bwa
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/nanopolish
    sbatch $ProgDir/sub_bwa_nanopolish.sh $Assembly $ReadDir/"$Strain"_concatenated_reads_filtered.fastq $OutDir/nanopolish
  done

  for Assembly in $(ls assembly/flye/F.venenatum/WT_minion/racon_10/WT_flye_racon10_renamed.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  echo "$Organism - $Strain"
  ReadDir=raw_dna/nanopolish/F.venenatum/WT_minion
  OutDir=$(dirname $Assembly)/nanopolish_bwa
  mkdir -p $OutDir
  ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/nanopolish
  sbatch $ProgDir/sub_bwa_nanopolish.sh $Assembly $ReadDir/"$Strain"_concatenated_reads_filtered.fastq $OutDir/nanopolish
done

  for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/racon_10/WT_miniasm_racon10_renamed.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  echo "$Organism - $Strain"
  ReadDir=raw_dna/nanopolish/F.venenatum/WT_minion
  OutDir=$(dirname $Assembly)/nanopolish_bwa
  mkdir -p $OutDir
  ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/nanopolish
  sbatch $ProgDir/sub_bwa_nanopolish.sh $Assembly $ReadDir/"$Strain"_concatenated_reads_filtered.fastq $OutDir/nanopolish
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_minion/racon_10/WT_minion_racon10_renamed.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  echo "$Organism - $Strain"
  OutDir=$(dirname $Assembly)/nanopolish
  RawReads=$(ls raw_dna/nanopolish/$Organism/$Strain/"$Strain"_concatenated_reads_filtered.fastq)
  AlignedReads=$(ls $OutDir/reads.sorted.bam)
  Ploidy=1
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
  sbatch $ProgDir/nanopolish_variants2.sh $Assembly $RawReads $AlignedReads $Ploidy $OutDir/variants
done
```

```bash
for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/racon_10/WT_miniasm_racon10_renamed.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  echo "$Organism - $Strain"
  OutDir=$(dirname $Assembly)/nanopolish_bwa/nanopolish
  RawReads=$(ls raw_dna/nanopolish/$Organism/$Strain/"$Strain"_concatenated_reads_filtered.fastq)
  AlignedReads=$(ls $OutDir/reads.sorted.bam)
  Ploidy=1
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
  sbatch $ProgDir/nanopolish_variants2.sh $Assembly $RawReads $AlignedReads $Ploidy $OutDir/variants
done

for Assembly in $(ls assembly/flye/F.venenatum/WT_minion/racon_10/WT_flye_racon10_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)/nanopolish_bwa/nanopolish
RawReads=$(ls raw_dna/nanopolish/$Organism/$Strain/"$Strain"_concatenated_reads_filtered.fastq)
AlignedReads=$(ls $OutDir/reads.sorted.bam)
Ploidy=1
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
sbatch $ProgDir/nanopolish_variants2.sh $Assembly $RawReads $AlignedReads $Ploidy $OutDir/variants
done
```








```bash
for Assembly in $(ls assembly/flye/F.venenatum/WT_minion/racon_10/WT_flye_racon10_renamed.fasta); do
Strain=WT_minion
Organism=F.venenatum
echo "$Organism - $Strain"
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
BuscoDB=$(ls -d /projects/oldhome/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd9
sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
done
```


```bash
ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
touch tmp.txt
for Assembly in $(ls assembly/flye/F.venenatum/WT_minion/racon_10/racon_round_10.fasta); do
OutDir=$(dirname $Assembly)
$ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/WT_flye_racon10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
rm tmp.txt
```
```bash
for Assembly in $(ls assembly/flye/F.venenatum/WT_minion/racon_10/WT_flye_racon10_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev | sed 's/_minion//g')
IlluminaDir=$(ls -d ../oldhome/groups/harrisonlab/project_files/fusarium_venenatum/qc_dna/paired/$Organism/$Strain)
echo $Strain
echo $Organism
TrimF2_Read=$(ls $IlluminaDir/F/FvenWT_S2_L001_R1_001_trim.fq.gz | head -n2 | tail -n1);
TrimR2_Read=$(ls $IlluminaDir/R/FvenWT_S2_L001_R2_001_trim.fq.gz | head -n2 | tail -n1);
TrimF3_Read=$(ls $IlluminaDir/F/FvenWT_S3_L001_R1_001_trim.fq.gz | head -n3 | tail -n1);
TrimR3_Read=$(ls $IlluminaDir/R/FvenWT_S3_L001_R2_001_trim.fq.gz | head -n3 | tail -n1);
echo $TrimF2_Read
echo $TrimR2_Read
echo $TrimF3_Read
echo $TrimR3_Read
OutDir=$(dirname $Assembly)
Iterations=10
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/pilon
sbatch $ProgDir/sub_pilon_2_libs.sh $Assembly $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir $Iterations
done
```
