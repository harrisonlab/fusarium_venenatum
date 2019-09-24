# Fv vs Fg

## Read alignment

Read alignments were performed to identify LS regions

### Downloading data

Fg illumina reads were downloaded fnrom Kelly and Ward 2018

NA1 population
37425
```bash
Organism="F.graminearum"
Strain="37425"
OutDir=/data/scratch/armita/fusarium/venenatum/raw_dna/paired/$Organism/$Strain
# fastq-dump --outdir $OutDir --split-files SRR5948915
fastq-dump --outdir $OutDir --split-3 SRR5948915
gzip $OutDir/*_1.fastq $OutDir/*_2.fastq
mkdir $OutDir/F
mkdir $OutDir/R
mv $OutDir/*_1.fastq.gz $OutDir/F/.
mv $OutDir/*_2.fastq.gz $OutDir/R/.
rm $OutDir/*.fastq
```

NA2 population

F335
```bash
Organism="F.graminearum"
Strain="F335"
OutDir=/data/scratch/armita/fusarium/venenatum/raw_dna/paired/$Organism/$Strain
fastq-dump --outdir $OutDir --split-3 SRR5948910
gzip $OutDir/*_1.fastq $OutDir/*_2.fastq
mkdir $OutDir/F
mkdir $OutDir/R
mv $OutDir/*_1.fastq.gz $OutDir/F/.
mv $OutDir/*_2.fastq.gz $OutDir/R/.
rm $OutDir/*.fastq
```


NA3 population

66031
```bash
Organism="F.graminearum"
Strain="66031"
CurDir=$PWD
OutDir=/data/scratch/armita/fusarium/venenatum/raw_dna/paired/$Organism/$Strain
fastq-dump --outdir $OutDir --split-3 SRR5948966
gzip $OutDir/*_1.fastq $OutDir/*_2.fastq
mkdir $OutDir/F
mkdir $OutDir/R
mv $OutDir/*_1.fastq.gz $OutDir/F/.
mv $OutDir/*_2.fastq.gz $OutDir/R/.
rm $OutDir/*.fastq
```

F. gerlachii

```bash
Organism="F.gerlachii"
Strain="38380"
OutDir=/data/scratch/armita/fusarium/venenatum/raw_dna/paired/$Organism/$Strain
mkdir -p $OutDir
fastq-dump --outdir $OutDir --split-3 SRR5948927
gzip $OutDir/*_1.fastq $OutDir/*_2.fastq
mkdir $OutDir/F
mkdir $OutDir/R
mv $OutDir/*_1.fastq.gz $OutDir/F/.
mv $OutDir/*_2.fastq.gz $OutDir/R/.
rm $OutDir/*.fastq
```

F. lousianense


```bash
Organism="F.louisianense"
Strain="54197"
OutDir=/data/scratch/armita/fusarium/venenatum/raw_dna/paired/$Organism/$Strain
mkdir $OutDir
fastq-dump --outdir $OutDir --split-3 SRR5948925
gzip $OutDir/*_1.fastq $OutDir/*_2.fastq
mkdir $OutDir/F
mkdir $OutDir/R
mv $OutDir/*_1.fastq.gz $OutDir/F/.
mv $OutDir/*_2.fastq.gz $OutDir/R/.
rm $OutDir/*.fastq
```

F. boothii

```bash
Organism="F.boothii"
Strain="F562"
OutDir=/data/scratch/armita/fusarium/venenatum/raw_dna/paired/$Organism/$Strain
fastq-dump --outdir $OutDir --split-3 SRR5948912
gzip $OutDir/*_1.fastq $OutDir/*_2.fastq
mkdir $OutDir/F
mkdir $OutDir/R
mv $OutDir/*_1.fastq.gz $OutDir/F/.
mv $OutDir/*_2.fastq.gz $OutDir/R/.
rm $OutDir/*.fastq
```



```bash
cd /data/scratch/armita/fusarium/venenatum
for StrainPath in $(ls -d raw_dna/paired/*/* | grep -v 'venenatum' | grep 'F335'); do
# Jobs=$(qstat | grep 'rna_qc_' | grep 'qw' | wc -l)
# while [ $Jobs -gt 1 ]; do
# sleep 1m
# printf "."
# Jobs=$(qstat | grep 'rna_qc_' | grep 'qw' | wc -l)
# done		
echo $StrainPath
Read_F=$(ls $StrainPath/F/*.fastq.gz)
Read_R=$(ls $StrainPath/R/*.fastq.gz)
IluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
done
```
<!--
Read trimming didnt work with fastqmcf as the reads werent all paired. As such, Trimmomatic was run.

```bash
screen -a
qlogin -l h=blacklace11
cd /data/scratch/armita/fusarium/venenatum
for StrainPath in $(ls -d raw_dna/paired/*/* | grep -v 'venenatum'); do
# Jobs=$(qstat | grep 'rna_qc_' | grep 'qw' | wc -l)
# while [ $Jobs -gt 1 ]; do
# sleep 1m
# printf "."
# Jobs=$(qstat | grep 'rna_qc_' | grep 'qw' | wc -l)
# done		
echo $StrainPath
Organism=$(echo $StrainPath | cut -f3 -d '/')
Strain=$(echo $StrainPath | cut -f4 -d '/')
Read_F=$(ls $StrainPath/F/*.fastq.gz)
Read_R=$(ls $StrainPath/R/*.fastq.gz)
IluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
OutDir=qc_dna/$Organism/$Strain
mkdir -p $OutDir/F
mkdir -p $OutDir/R
mkdir -p $OutDir/U
java -jar $HOME/prog/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 12 $Read_F $Read_R  $OutDir/F/${Strain}_trim_F.fastq.gz $OutDir/R/${Strain}_trim_R.fastq.gz $OutDir/U/${Strain}_trim_U1.fastq.gz $OutDir/U/${Strain}_trim_U2.fastq.gz ILLUMINACLIP:$IluminaAdapters:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
done
```
-->


### Fg vs Fv

LS regions were identified in Fv

```bash
cd /data/scratch/armita/fusarium/venenatum
for Reference in $(ls ../../../../../home/groups/harrisonlab/project_files/fusarium_venenatum/repeat_masked/F.venenatum/WT/illumina_assembly_ncbi/WT_contigs_unmasked.fa ../../../../../home/groups/harrisonlab/project_files/fusarium_venenatum/repeat_masked/F.venenatum/WT_minion/minion_submission/*_contigs_unmasked.fa | grep 'minion'); do
RefStrain=$(echo $Reference | rev | cut -f3 -d '/' | rev)
for StrainPath in $(ls -d raw_dna/paired/*/*); do
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
# F_Read=$(ls $StrainPath/F/*.fastq.gz)
# R_Read=$(ls $StrainPath/R/*.fastq.gz)
F_Read=$(ls $StrainPath/F/*.fastq.gz)
R_Read=$(ls $StrainPath/R/*.fastq.gz)
echo $F_Read
echo $R_Read
Prefix="${Organism}_${Strain}"
OutDir=analysis/genome_alignment/bwa/$Organism/$Strain/vs_${RefStrain}
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
qsub $ProgDir/sub_bwa.sh $Prefix $Reference $F_Read $R_Read $OutDir
done
done
```

```bash
cd /data/scratch/armita/fusarium/venenatum
for Reference in $(ls ../../../../../home/groups/harrisonlab/project_files/fusarium_venenatum/repeat_masked/F.venenatum/WT/illumina_assembly_ncbi/WT_contigs_unmasked.fa ../../../../../home/groups/harrisonlab/project_files/fusarium_venenatum/repeat_masked/F.venenatum/WT_minion/minion_submission/*_contigs_unmasked.fa | grep 'minion'); do
RefStrain=$(echo $Reference | rev | cut -f3 -d '/' | rev)
for StrainPath in $(ls -d ../../../../../home/groups/harrisonlab/project_files/fusarium_venenatum/qc_dna/paired/F.venenatum/WT); do
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz | tail -n1)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz | tail -n1)
echo $F_Read
echo $R_Read
Prefix="${Organism}_${Strain}"
OutDir=analysis/genome_alignment/bwa/$Organism/$Strain/vs_${RefStrain}
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
qsub $ProgDir/sub_bwa.sh $Prefix $Reference $F_Read $R_Read $OutDir
done
done
```

### Read coverage

Identify read coverage over each bp

```bash
  for Bam in $(ls analysis/genome_alignment/bwa/*/*/vs_*/*_sorted.bam | grep 'F335'); do
    Target=$(echo $Bam | rev | cut -f2 -d '/' | rev)
    Strain=$(echo $Bam | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Bam | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain - $Target"
    OutDir=$(dirname $Bam)
    samtools depth -aa $Bam > $OutDir/${Organism}_${Strain}_${Target}_depth.tsv
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_${Target}_depth.tsv > $OutDir/${Organism}_${Strain}_${Target}_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_${Target}_depth_10kb.tsv
  done
  for Target in "vs_WT" "vs_WT_minion"; do
    OutDir=analysis/genome_alignment/bwa/grouped_${Target}
    mkdir -p $OutDir
    cat analysis/genome_alignment/bwa/*/*/*/*_*_${Target}_depth_10kb.tsv > $OutDir/${Target}_grouped_depth.tsv
  done
```

### Plot read coverage


```R
library(readr)
setwd("~/Downloads/Aalt/coverage2")

appended_df <- read_delim("~/Downloads/Aalt/coverage2/vs_1166_grouped_depth.tsv", "\t", escape_double = FALSE, col_names = FALSE, col_types = cols(X4 = col_factor(levels = c("675", "97.0013", "97.0016", "650", "648", "24350", "1082", "1164", "635", "743", "1166", "1177"))), trim_ws = TRUE)

myFun <- function(x) {
  c(min = min(x), max = max(x),
    mean = mean(x), median = median(x),
    std = sd(x))
}

colnames(appended_df) <- c("contig","position", "depth", "strain")

appended_df$treatment <- paste(appended_df$strain , appended_df$contig)
tapply(appended_df$depth, appended_df$treatment, myFun)

df2 <- cbind(do.call(rbind, tapply(appended_df$depth, appended_df$treatment, myFun)))
write.csv(df2, '1166_contig_coverage.csv')

appended_df$depth <- ifelse(appended_df$depth > 100, 100, appended_df$depth)

# install.packages("ggplot2")
library(ggplot2)
require(scales)

for (i in 1:22){
contig = paste("contig", i, sep = "_")
p0 <- ggplot(data=appended_df[appended_df$contig == contig, ], aes(x=`position`, y=`depth`, group=1)) +
geom_line() +
labs(x = "Position (bp)", y = "Coverage") +
scale_y_continuous(breaks=seq(0,100,25), limits=c(0,100)) +
facet_wrap(~strain, nrow = 12, ncol = 1, strip.position = "left")
outfile = paste("1166_contig", i, "cov.jpg", sep = "_")
ggsave(outfile , plot = p0, device = 'jpg', path = NULL,
scale = 1, width = 500, height = 500, units = 'mm',
dpi = 150, limitsize = TRUE)
}

```

### Fv vs Fg

LS regions were identified in Fg


```bash


```


## Whole genome alignment

Nucmer was used to align Fv to Fg and Fg vs Fv

```bash
cd /home/groups/harrisonlab/project_files/fusarium_venenatum
Reference=$(ls assembly/external_group/F.graminearum/PH1/dna/Fusarium_graminearum.RR1.dna.toplevel.fa)
RefStrain=$(echo $Reference | rev | cut -f3 -d '/' | rev)
for Query in $(ls repeat_masked/F.venenatum/WT/illumina_assembly_ncbi/WT_contigs_unmasked.fa repeat_masked/F.venenatum/WT_minion/minion_submission/*_contigs_unmasked.fa); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix=${Strain}_vs_${RefStrain}
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```

```bash
cd /home/groups/harrisonlab/project_files/fusarium_venenatum
for Reference in $(ls repeat_masked/F.venenatum/WT/illumina_assembly_ncbi/WT_contigs_unmasked.fa repeat_masked/F.venenatum/WT_minion/minion_submission/*_contigs_unmasked.fa); do
RefStrain=$(echo $Reference | rev | cut -f3 -d '/' | rev)
for Query in $(ls assembly/external_group/F.graminearum/PH1/dna/Fusarium_graminearum.RR1.dna.toplevel.fa); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix=${Strain}_vs_${RefStrain}
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```
