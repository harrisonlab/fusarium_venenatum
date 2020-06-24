# Gene prediction minion/medaka genome assembly

## RepeatMasker

RepeatMasker screens DNA sequences for interspersed repeats and low complexity sequences.

### RepeatMasker and transposonPSI

```bash
conda activate general_tools

  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Repeat_masking
  BestAssembly=assembly/SMARTdenovo/F.venenatum/WT_minion/racon_10/medaka/medaka/pilon/WT_SMARTdenovo_medaka_pilon10_renamed.fasta
  OutDir=repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka
  sbatch $ProgDir/rep_modeling.sh $BestAssembly $OutDir
  sbatch $ProgDir/transposonPSI.sh $BestAssembly $OutDir
```

The TransposonPSI masked bases were used to mask additional bases from the repeatmasker / repeatmodeller softmasked and hardmasked files.


```bash
  for File in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked.fa); do
    OutDir=$(dirname $File)
    TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
    OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
    echo "$OutFile"
    bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
    echo "Number of masked bases:"
    cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
  done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
```

## Star

Spliced Transcripts Alignment to a Reference. 


```bash
conda install star

#samtools sort used for sorting
#for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_unmasked.fa)
#do
#Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev) 
#Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
#echo $Assembly
#echo "$Organism - $Strain"
#for FileF in $(ls -d ../oldhome/groups/harrisonlab/project_files/quorn/filtered/WTCHG_*.1.fq)
#do
#FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/1.fq/2.fq/g')
#echo $FileF
#echo $FileR
#Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/.1.fq//g')
#OutDir=alignment/star/$Organism/$Strain/medaka_assembly_samtoolsort/$Sample_Name
#ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
#sbatch $ProgDir/star.sh $Assembly $FileF $FileR $OutDir
#done
#done

#Star sortbycoord used
  for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev) 
    Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
    echo $Assembly
    echo "$Organism - $Strain"
    for FileF in $(ls -d ../oldhome/groups/harrisonlab/project_files/quorn/filtered/WTCHG_258647_224.1.fq) ; do
      FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/1.fq/2.fq/g')
      echo $FileF
      echo $FileR
      #Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
      #echo "$Timepoint"
      Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/.1.fq//g')
      OutDir=alignment/star/$Organism/$Strain/medaka_assembly_star/$Sample_Name
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
      sbatch $ProgDir/star2.sh $Assembly $FileF $FileR $OutDir
    done
  done
```

If multiple RNAseq samples are used, alignment outputs can be concatenated using samtools. 

```bash
# For all alignments
#BamFiles=$(ls alignment/star/F.venenatum/WT_minion/medaka_assembly/*/*.sorted.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
#OutDir=/data/scratch/gomeza/fusarium_venenatum/star/F.venenatum/WT_minion/medaka_assembly/concatenated
#mkdir -p $OutDir
#samtools merge -f $OutDir/concatenated.bam $BamFiles

#BamFiles=$(ls ../oldhome/groups/harrisonlab/project_files/fusarium_venenatum/alignment/star/F.venenatum/WT/treatment/*/*.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
#OutDir=/data/scratch/gomeza/fusarium_venenatum/previous_alignment/star/F.venenatum/WT/concatenated
#mkdir -p $OutDir
#samtools merge -f $OutDir/concatenated.bam $BamFiles

# One block
#BamFiles=$(ls alignment/star/F.venenatum/WT_minion/medaka_assembly/Rep1/*/*.sorted.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
#OutDir=alignment/star/F.venenatum/WT_minion/medaka_assembly/concatenated
#mkdir -p $OutDir
#samtools merge -f $OutDir/concatenated_258646_1rep.bam $BamFiles

# One group
#BamFiles=$(ls alignment/star/F.venenatum/WT_minion/medaka_assembly/WTCHG_258647*/*.sorted.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
#OutDir=alignment/star/F.venenatum/WT_minion/medaka_assembly/concatenated_258647
#mkdir -p $OutDir
#samtools merge -f $OutDir/concatenated_258647.bam $BamFiles

#BamFiles=$(ls alignment/star/F.venenatum/WT_minion/medaka_assembly_star/*/*.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
#OutDir=alignment/star/F.venenatum/WT_minion/medaka_assembly_v2/concatenated
#mkdir -p $OutDir
#samtools merge -f $OutDir/concatenated_259732.bam $BamFiles

# WT media samples only. This will be used for Codingquarry
BamFiles=$(ls alignment/star/F.venenatum/WT_minion/medaka_assembly_star/one_per_media/*/*.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
OutDir=alignment/star/F.venenatum/WT_minion/medaka_assembly_star/one_per_media/concatenated
mkdir -p $OutDir
samtools merge -f $OutDir/concatenated.bam $BamFiles

#All alignments. This will be used for Braker
screen -a
BamFiles=$(ls alignment/star/F.venenatum/WT_minion/medaka_assembly_star/*/*.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
OutDir=/data/scratch/gomeza/fusarium_venenatum/star/F.venenatum/WT_minion/concatenated
mkdir -p $OutDir
samtools merge -f $OutDir/concatenated.bam $BamFiles
```

## Braker

Braker is a combination of GeneMark-ET and Augustus used for gene annotation.

```bash
conda activate gene_pred

#for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
#Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev) 
#Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
#echo "$Organism - $Strain"
#OutDir=gene_pred/braker/$Organism/$Strain
#AcceptedHits=../../data/scratch/gomeza/fusarium_venenatum/star/F.venenatum/WT_minion/medaka_assembly/concatenated/concatenated.bam
#GeneModelName="$Organism"_"$Strain"_braker_all
#ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
#sbatch $ProgDir/braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
#done

#for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
#Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev) 
#Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
#echo "$Organism - $Strain"
#OutDir=gene_pred/braker/$Organism/$Strain/Rep1
#AcceptedHits=alignment/star/F.venenatum/WT_minion/medaka_assembly/concatenated/concatenated_258646_1rep.bam
#GeneModelName="$Organism"_"$Strain"_braker_medaka_rep1
#ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
#sbatch $ProgDir/braker_fungi_v2.sh $Assembly $OutDir $AcceptedHits $GeneModelName
#done

#for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
#Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev) 
#Organism=$(echo $Assembly | rev | cut -d '/' -f5 | rev)
#echo "$Organism - $Strain"
#OutDir=gene_pred/braker/$Organism/$Strain/group_47
#AcceptedHits=alignment/star/F.venenatum/WT_minion/medaka_assembly/concatenated_258647/concatenated_258647.bam
#GeneModelName="$Organism"_"$Strain"_braker_medaka_47_5
#ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
#sbatch $ProgDir/braker_fungi_v2.sh $Assembly $OutDir $AcceptedHits $GeneModelName
#done

#for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
#Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev) 
#Organism=$(echo $Assembly | rev | cut -d '/' -f5 | rev)
#echo "$Organism - $Strain"
#OutDir=gene_pred/braker_v2/$Organism/$Strain/
#AcceptedHits=alignment/star/F.venenatum/WT_minion/medaka_assembly_star/concatenated/concatenated_259732.bam
#GeneModelName="$Organism"_"$Strain"_braker_medaka_32
#ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
#sbatch $ProgDir/braker_fungi_v2.sh $Assembly $OutDir $AcceptedHits $GeneModelName
#done
```


## Genome-guided assembly and CodingQuarry

CodingQuarry in pathogen mode is used to predict aditional genes and added to braker predictions


#### Stringtie RNA-seq alignments assembler

```bash
#for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_unmasked.fa); do
#Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev) 
#Organism=$(echo $Assembly| rev | cut -d '/' -f5 | rev) 
#echo "$Organism - $Strain"
#OutDir=gene_pred/stringtie/$Organism/$Strain/concatenated_prelim_46_1rep
#mkdir -p $OutDir
#AcceptedHits=alignment/star/F.venenatum/WT_minion/medaka_assembly/concatenated/concatenated_258646_1rep.bam
#ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
#sbatch $ProgDir/stringtie.sh $AcceptedHits $OutDir
#done

#for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_unmasked.fa); do
#Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev) 
#Organism=$(echo $Assembly| rev | cut -d '/' -f5 | rev) 
#echo "$Organism - $Strain"
#OutDir=gene_pred/stringtie_v2/$Organism/$Strain/concatenated_prelim_32
#mkdir -p $OutDir
#AcceptedHits=alignment/star/F.venenatum/WT_minion/medaka_assembly_star/concatenated/concatenated_259732.bam
#ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
#sbatch $ProgDir/stringtie.sh $AcceptedHits $OutDir
#done
```

#### Cufflinks RNA-seq alignments assembler

```bash
#for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_unmasked.fa); do
#Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev) 
#Organism=$(echo $Assembly| rev | cut -d '/' -f5 | rev) 
#echo "$Organism - $Strain"
#OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim_47
#mkdir -p $OutDir
#AcceptedHits=alignment/star/F.venenatum/WT_minion/medaka_assembly/concatenated_258647/concatenated_258647.bam
#ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
#sbatch $ProgDir/cufflinks.sh $AcceptedHits $OutDir
#done

#for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_unmasked.fa); do
#Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev) 
#Organism=$(echo $Assembly| rev | cut -d '/' -f5 | rev) 
#echo "$Organism - $Strain"
#OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim_all
#mkdir -p $OutDir
#AcceptedHits=../../data/scratch/gomeza/fusarium_venenatum/star/F.venenatum/WT_minion/medaka_assembly/concatenated/concatenated.bam
#ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
#sbatch $ProgDir/cufflinks2.sh $AcceptedHits $OutDir
#done

#or Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_unmasked.fa); do
#Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev) 
#Organism=$(echo $Assembly| rev | cut -d '/' -f5 | rev) 
#echo "$Organism - $Strain"
#OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim_Rep1
#mkdir -p $OutDir
#AcceptedHits=alignment/star/F.venenatum/WT_minion/medaka_assembly/concatenated/concatenated_258646_1rep.bam
#ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
#sbatch $ProgDir/cufflinks.sh $AcceptedHits $OutDir
#done


#for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_unmasked.fa); do
#Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev) 
#Organism=$(echo $Assembly| rev | cut -d '/' -f5 | rev) 
#echo "$Organism - $Strain"
#OutDir=gene_pred/cufflinks_v2/$Organism/$Strain/concatenated_prelim_32
#mkdir -p $OutDir
#AcceptedHits=alignment/star/F.venenatum/WT_minion/medaka_assembly_v2/concatenated/concatenated_259732.bam
#ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
#sbatch $ProgDir/cufflinks.sh $AcceptedHits $OutDir
#done
```

## Cufflinks RNA-seq alignments assembler

```bash
  for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev) 
    Organism=$(echo $Assembly| rev | cut -d '/' -f5 | rev) 
    echo "$Organism - $Strain"
    OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim_vAG
    mkdir -p $OutDir
    AcceptedHits=alignment/star/F.venenatum/WT_minion/one_per_media/concatenated/concatenated.bam
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
    sbatch $ProgDir/cufflinks.sh $AcceptedHits $OutDir
  done
```

## CodinQuarry 

Note: run_CQ-PM_stranded.sh and run_CQ-PM_unstranded.sh scripts are included in cndigquarry scripts are used to run CQ pathogen mode using signalp 4.1. The script in this folder was edited to use signalp5. 

```bash
  #for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_unmasked.fa); do
    #Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev) 
    #Organism=$(echo $Assembly| rev | cut -d '/' -f5 | rev)
    #echo "$Organism - $Strain"
    #OutDir=gene_pred/codingquary/$Organism/$Strain/Rep1
    #mkdir -p $OutDir
    #GTF=gene_pred/stringtie/F.venenatum/WT_minion/concatenated_prelim_Rep1/out.gtf
    #ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
    #sbatch $ProgDir/codingquarry.sh $Assembly $GTF $OutDir
  #done

  #for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_unmasked.fa); do
  #Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev) 
  #Organism=$(echo $Assembly| rev | cut -d '/' -f5 | rev)
  #echo "$Organism - $Strain"
  #OutDir=gene_pred/codingquary/$Organism/$Strain/Rep1_sig4
  #mkdir -p $OutDir
  #GTF=gene_pred/stringtie/F.venenatum/WT_minion/concatenated_prelim_Rep1/out.gtf
  #ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
  #sbatch $ProgDir/codingquarry2.sh $Assembly $GTF $OutDir
#done
```
```bash
#for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_unmasked.fa); do
  #Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev) 
  #Organism=$(echo $Assembly| rev | cut -d '/' -f5 | rev) 
  #echo "$Organism - $Strain"
  #OutDir=gene_pred/codingquary_cuff/$Organism/$Strain/47
  #mkdir -p $OutDir
  #GTF=gene_pred/cufflinks/F.venenatum/WT_minion/concatenated_prelim_47/transcripts.gtf
  #ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
  #sbatch $ProgDir/codingquarry2.sh $Assembly $GTF $OutDir
#done

#for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_unmasked.fa); do
  #Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev) 
  #Organism=$(echo $Assembly| rev | cut -d '/' -f5 | rev) 
  #echo "$Organism - $Strain"
  #OutDir=gene_pred/codingquary_cuff/$Organism/$Strain/Rep1
  #mkdir -p $OutDir
  #GTF=gene_pred/cufflinks/F.venenatum/WT_minion/concatenated_prelim_Rep1/transcripts.gtf
  #ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
  #sbatch $ProgDir/codingquarry2.sh $Assembly $GTF $OutDir
#done
```

## Codingquarry 

```bash
conda activate antismash_py27

for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_unmasked.fa); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev) 
  Organism=$(echo $Assembly| rev | cut -d '/' -f5 | rev) 
  echo "$Organism - $Strain"
  OutDir=gene_pred/codingquarry/$Organism/$Strain
  mkdir -p $OutDir
  GTF=gene_pred/cufflinks/F.venenatum/WT_minion/concatenated_prelim_vAG/transcripts.gtf
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
  sbatch $ProgDir/codingquarry2.sh $Assembly $GTF $OutDir
done
```

#### Add additional transcripts to Braker gene models.


Additional transcripts predicted by CodingQuarry are added to the final gene models.

```bash
  # The following perl scripts requires the installation of some libraries. Run these commands in a perly environment.
  # Install the required libraries (if any) using cpanm
  # cpanm Bio::Perl

BrakerGff=gene_pred/braker/F.venenatum/WT_minion/group_47/augustus.hints.gff3
Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev)
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
Assembly=$(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
CodingQuarryGff=gene_pred/codingquary_cuff/F.venenatum/WT_minion/47/out/PredictedPass.gff3
PGNGff=gene_pred/codingquary_cuff/F.venenatum/WT_minion/47/out/PGN_predictedPass.gff3
AddDir=gene_pred/codingquary_cuff/$Organism/$Strain/47/additional
FinalDir=gene_pred/codingquary_cuff/$Organism/$Strain/47/final
AddGenesList=$AddDir/additional_genes.txt
AddGenesGff=$AddDir/additional_genes.gff
FinalGff=$AddDir/combined_genes.gff
mkdir -p $AddDir
mkdir -p $FinalDir

  # Create a list with the additional transcripts in CondingQuarry gff (and CQPM) vs Braker gene models
bedtools intersect -v -a $CodingQuarryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
  
  # Creat Gff file with the additional transcripts
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuarryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
	
  # Create a final Gff file with gene features
$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3

# Create fasta files from each gene feature in the CodingQuarry gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary

# Create fasta files from each gene feature in the Braker gff3
cp $BrakerGff $FinalDir/final_genes_Braker.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker

# Combine both fasta files
cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

# Combine both gff3 files
GffBraker=$FinalDir/final_genes_CodingQuary.gff3
GffQuary=$FinalDir/final_genes_Braker.gff3
GffAppended=$FinalDir/final_genes_appended.gff3
cat $GffBraker $GffQuary > $GffAppended

# Check the final number of genes

for DirPath in $(ls -d $FinalDir); do
echo $DirPath;
cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
echo "";
done
```

```
gene_pred/codingquary_cuff/F.venenatum/WT_minion/47/final
12529
1415
13944
```


  #### Remove duplicate and rename genes.

  ```bash
  for GffAppended in $(ls gene_pred/codingquary_cuff/F.venenatum/WT_minion/47/final/final_genes_appended.gff3);
  do
    Strain=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
    Organism=$(echo $GffAppended | rev | cut -d '/' -f5 | rev)
    echo "$Organism - $Strain"
    FinalDir=gene_pred/codingquary_cuff/F.venenatum/WT_minion/47/final
    # Remove duplicated genes
    GffFiltered=gene_pred/codingquary_cuff/F.venenatum/WT_minion/47/final/filtered_duplicates.gff
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
    $ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
    # Rename genes
    GffRenamed=gene_pred/codingquary_cuff/F.venenatum/WT_minion/47/final/final_genes_appended_renamed.gff3
    LogFile=gene_pred/codingquary_cuff/F.venenatum/WT_minion/47/final/final_genes_appended_renamed.log
    $ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
    rm $GffFiltered
    # Create renamed fasta files from each gene feature   
    Assembly=$(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
    $ProgDir/gff2fasta.pl $Assembly $GffRenamed gene_pred/codingquary_cuff/F.venenatum/WT_minion/47/final/final_genes_appended_renamed
    # The proteins fasta file contains * instead of Xs for stop codons, these should be changed
    sed -i 's/\*/X/g' gene_pred/codingquary_cuff/F.venenatum/WT_minion/47/final/final_genes_appended_renamed.pep.fasta
  done 
```

## Busco

```bash
  for Assembly in $(ls gene_pred/codingquary_cuff/F.venenatum/WT_minion/47/final/final_genes_appended_renamed.gene.fasta); do
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
    BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
    OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10
    sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
  done
```

## Interproscan

Interproscan was used to give gene models functional annotations.


```bash
# This command will split your gene fasta file and run multiple interproscan jobs.
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  for Genes in $(ls gene_pred/codingquary_cuff/F.venenatum/WT_minion/47/final/final_genes_appended_renamed.pep.fasta); do
    echo $Genes
    $ProgDir/interproscan.sh $Genes
  done 2>&1 | tee -a interproscan_submisison.log
```




Following interproscan annotation split files were combined using the following commands:

```bash
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  for Proteins in $(ls gene_pred/codingquarry_cuff_final/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta); do
    Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    echo $Strain
    InterProRaw=gene_pred/interproscan/AG_medaka_47/raw
    $ProgDir/append_interpro.sh $Proteins $InterProRaw
  done
```

## B) SwissProt

```bash
for Proteome in $(ls gene_pred/codingquarry_cuff_final/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/swissprot/$Organism/$Strain
SwissDbDir=../dbUniprot/swissprot_2020_June
SwissDbName=uniprot_sprot.db
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
sbatch $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
done
```
```

## Effector genes

Putative pathogenicity and effector related genes were identified within Braker
gene models using a number of approaches:

 * A) From Augustus gene models - Identifying secreted proteins
 * B) From Augustus gene models - Effector identification using EffectorP


### A) From Augustus gene models - Identifying secreted proteins

 Required programs:
  * SignalP-4.1
  * TMHMM

 Proteins that were predicted to contain signal peptides were identified using
 the following commands:


```bash
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
CurPath=$PWD
for Proteome in $(ls gene_pred/codingquarry_cuff_final/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
SplitDir=gene_pred/final_genes_split/$Organism/$Strain
mkdir -p $SplitDir
BaseName="$Organism""_$Strain"_final_preds
$ProgDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
for File in $(ls $SplitDir/*_final_preds_*); do
sbatch $ProgDir/pred_signalP.sh $File signalp-4.1
done
done
```

 The batch files of predicted secreted proteins needed to be combined into a
 single file for each strain. This was done with the following commands:
 ```bash
for SplitDir in $(ls -d gene_pred/final_genes_split/F.venenatum/WT_minion); do
Strain=$(echo $SplitDir | rev |cut -d '/' -f1 | rev)
Organism=$(echo $SplitDir | rev |cut -d '/' -f2 | rev)
InStringAA=''
InStringNeg=''
InStringTab=''
InStringTxt=''
SigpDir=final_genes_signalp-4.1
for GRP in $(ls -l $SplitDir/*_final_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do
InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.aa";
InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp_neg.aa";
InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.tab";
InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.txt";
done
cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.aa
cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_neg_sp.aa
tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.tab
cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.txt
done
done
 ```


 Some proteins that are incorporated into the cell membrane require secretion.
 Therefore proteins with a transmembrane domain are not likely to represent
 cytoplasmic or apoplastic effectors.


PATH=${PATH}:/data/scratch/gomeza/prog/tmhmm-2.0c/bin


 Proteins containing a transmembrane domain were identified:

 ```bash
for Proteome in $(ls gene_pred/codingquarry_cuff_final/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
sbatch $ProgDir/TMHMM.sh $Proteome
done
 ```

 Those proteins with transmembrane domains were removed from lists of Signal peptide containing proteins

 ```bash
 for File in $(ls gene_pred/trans_mem/F.venenatum/WT_minion/WT_minion_TM_genes_neg.txt); do
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  echo "$Organism - $Strain"
  TmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
  cat $File | cut -f1 > $TmHeaders
  SigP=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp.aa)
  OutDir=$(dirname $SigP)
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  $ProgDir/extract_from_fasta.py --fasta $SigP --headers $TmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
  cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l
 done
```

## B) From Augustus gene models - Effector identification using EffectorP

Required programs:
 * EffectorP.py

```bash
for Proteome in $(ls gene_pred/codingquarry_cuff_final/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
BaseName="$Organism"_"$Strain"_EffectorP
OutDir=analysis/effectorP/$Organism/$Strain
EffectorP.py -o "$BaseName".txt -E "$BaseName".fa -i $Proteome
mv "$BaseName".txt $OutDir
mv "$BaseName".fa $OutDir
#ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
#sbatch $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
done
```

```bash
for File in $(ls analysis/effectorP/*/*/*_EffectorP.txt); do
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  echo "$Organism - $Strain"
  Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_headers.txt/g')
  cat $File | grep 'Effector' | cut -f1 > $Headers
  Secretome=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp_no_trans_mem.aa)
  OutFile=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.aa/g')
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  $ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > $OutFile
  OutFileHeaders=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted_headers.txt/g')
  cat $OutFile | grep '>' | tr -d '>' > $OutFileHeaders
  cat $OutFileHeaders | wc -l
  Gff=$(ls gene_pred/codingquarry_cuff_final/F.venenatum/WT_minion/final/final_genes_appended_renamed.gff3)
  EffectorP_Gff=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.gff/g')
  $ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
  cat $EffectorP_Gff | grep -w 'gene' | wc -l
done > tmp.txt
```


```bash
for AntiSmash in $(ls analysis/secondary_metabolites/antismash/F.venenatum/WT_minion_vAG/*appended.gbk); do
Organism=$(echo $AntiSmash | rev | cut -f3 -d '/' | rev)
Strain=$(echo $AntiSmash | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
Prefix=$OutDir/WT_antismash_results
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
$ProgDir/antismash2gffv5.py --inp_antismash $AntiSmash --out_prefix $Prefix

# Identify secondary metabolites within predicted clusters
printf "Number of secondary metabolite detected:\t"
cat "$Prefix"_secmet_clusters.gff | wc -l
GeneGff=gene_pred/codingquarry_cuff_final/F.venenatum/WT_minion/final/final_genes_appended_renamed.gff3
bedtools intersect -u -a $GeneGff -b "$Prefix"_secmet_clusters.gff > "$Prefix"_secmet_genes.gff
cat "$Prefix"_secmet_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > "$Prefix"_antismash_secmet_genes.txt
bedtools intersect -wo -a $GeneGff -b "$Prefix"_secmet_clusters.gff | grep 'mRNA' | cut -f9,10,12,18 | sed "s/ID=//g" | perl -p -e "s/;Parent=g\w+//g" | perl -p -e "s/;Notes=.*//g" > "$Prefix"_secmet_genes.tsv
printf "Number of predicted proteins in secondary metabolite clusters:\t"
cat "$Prefix"_secmet_genes.txt | wc -l
printf "Number of predicted genes in secondary metabolite clusters:\t"
cat "$Prefix"_secmet_genes.gff | grep -w 'gene' | wc -l
done


SMURF was also run to identify secondary metabolite gene clusters.

Genes needed to be parsed into a specific tsv format prior to submission on the SMURF webserver.

Gff=$(ls gene_pred/codingquarry_cuff_final/F.venenatum/WT_minion/final/final_genes_appended_renamed.gff3)
OutDir=analysis/secondary_metabolites/smurf/F.venenatum/WT_minion
mkdir -p $OutDir
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
$ProgDir/gff2smurf.py --gff $Gff > $OutDir/WT_minion_genes_smurf.tsv
SMURF output was received by email and downloaded to the cluster in the output directory above.

Output files were parsed into gff format:

  OutDir=analysis/secondary_metabolites/smurf/F.venenatum/WT
  Prefix="WT"
  GeneGff=gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.gff3
  SmurfClusters=$OutDir/Secondary-Metabolite-Clusters.txt
  SmurfBackbone=$OutDir/Backbone-genes.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
  $ProgDir/smurf2gff.py --smurf_clusters $SmurfClusters --smurf_backbone $SmurfBackbone > $OutDir/Smurf_clusters.gff
  bedtools intersect -wo -a $GeneGff -b $OutDir/Smurf_clusters.gff | grep 'mRNA' | cut -f9,10,12,18 | sed "s/ID=//g" | perl -p -i -e "s/;Parent=g\w+//g" | perl -p -i -e "s/;Notes=.*//g" > $OutDir/"$Prefix"_smurf_secmet_genes.tsv
Total number of secondary metabolite clusters:

for Assembly in $(ls repeat_masked/*/*/illumina_assembly_ncbi/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w 'WT'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
mkdir -p $OutDir
GeneGff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
AntismashClusters=$(ls analysis/secondary_metabolites/antismash/$Organism/$Strain/*_secmet_clusters.gff)
SmurfClusters=$(ls analysis/secondary_metabolites/smurf/$Organism/$Strain/Smurf_clusters.gff)
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

F) Genes with transcription factor annotations:
A list of PFAM domains, superfamily annotations used as part of the DBD database and a further set of interproscan annotations listed by Shelest et al 2017 were made http://www.transcriptionfactor.org/index.cgi?Domain+domain:all https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5415576/

  for Interpro in $(ls gene_pred/interproscan/AG_medaka_47/F.venenatum/WT_minion/WT_minion_interproscan.tsv); do
    Organism=$(echo $Interpro | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Interpro | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/transcription_factors/vAG/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
    $ProgDir/interpro2TFs.py --InterPro $Interpro > $OutDir/"$Strain"_TF_domains.tsv
    echo "total number of transcription factors"
    cat $OutDir/"$Strain"_TF_domains.tsv | cut -f1 | sort | uniq > $OutDir/"$Strain"_TF_gene_headers.txt
    cat $OutDir/"$Strain"_TF_gene_headers.txt | wc -l
  done




  faidx -d '|' final_genes_appended_renamed.cdna.fasta $(tr '\n' ' ' < Tri5_genes.txt ) > Tri5_genes.fasta

  