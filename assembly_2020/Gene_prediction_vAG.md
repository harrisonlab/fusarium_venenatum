# RepeatMasker

RepeatMasker screens DNA sequences for interspersed repeats and low complexity sequences.

### Requirements

```bash
conda activate general_tools
```

## RepeatMasker and transposonPSI


```bash
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Repeat_masking
BestAssembly=assembly/SMARTdenovo/F.venenatum/WT_minion/racon_10/medaka/medaka/pilon/WT_SMARTdenovo_medaka_pilon10_renamed.fasta
OutDir=repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka
sbatch $ProgDir/rep_modeling.sh $BestAssembly $OutDir
sbatch $ProgDir/transposonPSI.sh $BestAssembly $OutDir
```


The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and hardmasked files.


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

### Requirements

```bash
conda install star
```

### Typical run

```bash
for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_unmasked.fa)
do
Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev) 
Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
echo $Assembly
echo "$Organism - $Strain"
for FileF in $(ls -d ../oldhome/groups/harrisonlab/project_files/quorn/filtered/WTCHG_259732_224.1.fq)
do
FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/1.fq/2.fq/g')
echo $FileF
echo $FileR
#Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
#echo "$Timepoint"
Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/.1.fq//g')
OutDir=alignment/star/$Organism/$Strain/medaka_assembly/$Sample_Name
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
sbatch $ProgDir/star.sh $Assembly $FileF $FileR $OutDir
done
done
```

If multiple RNAseq samples are used, alignment outputs can be concatenated using samtools. 

```bash
# For all alignments
BamFiles=$(ls alignment/star/F.venenatum/WT_minion/medaka_assembly/*/*.sorted.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
OutDir=/data/scratch/gomeza/fusarium_venenatum/star/F.venenatum/WT_minion/medaka_assembly/concatenated
mkdir -p $OutDir
samtools merge -f $OutDir/concatenated.bam $BamFiles

BamFiles=$(ls ../oldhome/groups/harrisonlab/project_files/fusarium_venenatum/alignment/star/F.venenatum/WT/treatment/*/*.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
OutDir=/data/scratch/gomeza/fusarium_venenatum/previous_alignment/star/F.venenatum/WT/concatenated
mkdir -p $OutDir
samtools merge -f $OutDir/concatenated.bam $BamFiles
```

# Gene prediction 

## Braker

Braker is a combination of GeneMark-ET and Augustus used for gene annotation.

### Requirements

Braker requires the instalation of multiple perl libraries, Augustus and Genemark. There are conda packages for Augustus and Braker. However, these packages are often lagging behing the releasess of AUGUSTUS and BRAKER. Therefore manual installation might be necessary.

### Conda installation

```bash
conda activate gene_pred # e.g.

for Assembly in $(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev) 
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/braker/$Organism/$Strain
AcceptedHits=../../data/scratch/gomeza/fusarium_venenatum/star/F.venenatum/WT_minion/medaka_assembly/concatenated/concatenated.bam
GeneModelName="$Organism"_"$Strain"_braker_medaka
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
sbatch $ProgDir/braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done

    for Assembly in $(ls assembly/previous_versions/F.venenatum/WT/WT_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev) 
    Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/braker/$Organism/$Strain
    AcceptedHits=path/to/your/spliced/aligments/files/*_aligmentAligned.sortedByCoord.out.bam # STAR output, see Genome_aligners folder
    #AcceptedHits=alignment/concatenated.bam # Concatenatented alignment files can be used
    GeneModelName="$Organism"_"$Strain"_braker 
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
    sbatch $ProgDir/braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  done
```


## Genome-guided assembly and CodingQuarry

CodingQuarry in pathogen mode is used to predict aditional genes and added to braker predictions

### Requirements

```bash
# Conda installation

# CodingQuarry requires a conda environment with python 2.7
# e.g. conda create --name gene_pred_py27 python=2.7

conda install stringtie
conda install codingquarry

# The environmental variable QUARRY_PATH is set in your profile (needed for CodingQuarry)

nano ~/.profile
export QUARRY_PATH="/home/"USER_ID"/miniconda3/envs/"USER_ENV_py27"/opt/codingquarry-2.0/QuarryFiles/QuarryFiles" #Add

# SignalP is needed. Add this path to your profile or 
vPATH=${PATH}:/data/scratch/gomeza/prog/signalp/signalp-5.0b/bin

. ~/.profile # Refresh your profile
```

### Typical run


#### Stringtie RNA-seq alignments assembler

```bash
  for Assembly in $(ls assembly/previous_versions/F.venenatum/WT/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f5 | rev) # Edit to set your ouput directory
    Organism=$(echo $Assembly| rev | cut -d '/' -f6 | rev) # Edit to set your ouput directory
    echo "$Organism - $Strain"
    OutDir=gene_pred/stringtie/$Organism/$Strain/concatenated_prelim
    mkdir -p $OutDir
    AcceptedHits=path/to/your/spliced/aligments/files.bam
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
    qbatch $ProgDir/stringties.sh $AcceptedHits $OutDir
   done
```


#### CodinQuarry 

Note: run_CQ-PM_stranded.sh and run_CQ-PM_unstranded.sh scripts are included in cndigquarry scripts are used to run CQ pathogen mode using signalp 4.1. The script in this folder was edited to use signalp5. 

```bash
  for Assembly in $(ls path/to/unmasked/genome/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f5 | rev) # Edit to set your ouput directory
    Organism=$(echo $Assembly| rev | cut -d '/' -f6 | rev) # Edit to set your ouput directory
    echo "$Organism - $Strain"
    OutDir=gene_pred/codingquary/$Organism/$Strain/
    mkdir -p $OutDir
    GTF=path/to/RNAseq/alignment/assembly/*.gtf # GFT file from stringtie/cufflinks output. See Genome-guided_assemblers scripts
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
    sbatch $ProgDir/codingquarry.sh $Assembly $GTF $OutDir
  done
```

#### Add additional transcripts to Braker gene models.


Additional transcripts predicted by CodingQuarry are added to the final gene models.

```bash
  # The following perl scripts requires the installation of some libraries. Run these commands in a perly environment.
  # Install the required libraries (if any) using cpanm
  # cpanm Bio::Perl

  BrakerGff=$(ls path/to/braker/gene/models/augustus.hints.gff3)
	Strain=$(echo $BrakerGff| rev | cut -d '/' -f2 | rev)
	Organism=$(echo $BrakerGff | rev | cut -d '/' -f3 | rev)
	echo "$Organism - $Strain"
	Assembly=$(ls path/to/softmasked/genome/assembly/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
	CodingQuarryGff=path/to/codingquarry/gff3/$Organism/$Strain/out/PredictedPass.gff3
	PGNGff=path/to/codingquarry/pathogen/mode/gff3/$Organism/$Strain/out/PGN_predictedPass.gff3
	AddDir=gene_pred/codingquary/$Organism/$Strain/additional # Additional transcripts directory
	FinalDir=gene_pred/codingquary/$Organism/$Strain/final # Final directory
	AddGenesList=$AddDir/additional_genes.txt
	AddGenesGff=$AddDir/additional_genes.gff
	FinalGff=$AddDir/combined_genes.gff
	mkdir -p $AddDir
	mkdir -p $FinalDir

  # Create a list with the additional transcripts in CondingQuarry gff (and CQPM) vs Braker gene models
	bedtools intersect -v -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
	bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
  
  # Creat Gff file with the additional transcripts
	ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
	$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
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

  #### Remove duplicate and rename genes.

  ```bash
  for GffAppended in $(ls $FinalDir/final_genes_appended.gff3);
  do
    Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    FinalDir=path/to/final/gene/predction/folder #/final
    # Remove duplicated genes
    GffFiltered=$FinalDir/filtered_duplicates.gff
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
    $ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
    # Rename genes
    GffRenamed=$FinalDir/final_genes_appended_renamed.gff3
    LogFile=$FinalDir/final_genes_appended_renamed.log
    $ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
    rm $GffFiltered
    # Create renamed fasta files from each gene feature   
    Assembly=$(ls assembly_vAG/canu_1step/N.ditissima/R0905/polished/repeat_masked/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
    $ProgDir/gff2fasta.pl $Assembly $GffRenamed $FinalDir/final_genes_appended_renamed
    # The proteins fasta file contains * instead of Xs for stop codons, these should be changed
    sed -i 's/\*/X/g' $FinalDir/final/final_genes_appended_renamed.pep.fasta
  done 
```