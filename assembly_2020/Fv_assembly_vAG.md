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
  for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/WT_minion.fa); do
    Strain=WT
    Organism=F.venenatum
    echo "$Organism - $Strain"
    ProgDir=/home/gomeza/git_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
    OutDir=$(dirname $Assembly)/busco
    sbatch $ProgDir/sub_busco.sh $Assembly $BuscoDB $OutDir
  done
```



Kat tool was used to assess level of error and duplications in the genome assemblies generated. Mapleson et al., 2016.

```bash
#Editting....
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
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/*10.fasta); do
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/*10.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
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






```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.txt
  for Assembly in $(ls assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/*10.fasta); do
    OutDir=$(dirname $Assembly)
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/Hg199_racon10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
  done
  rm tmp.txt
```
