# RNA-Seq analysis 

## Perform qc on RNA-Seq data

```bash
# Run fastqc
    for RawData in $(ls ../../../archives/2021_camb_general/20210128_Fvenenatum_CarbonRNAseq/rawdata/*.fq.gz | grep 'FvC0'); do
        echo $RawData
        Timepoint=$(echo $RawData | rev | cut -d '/' -f1 | rev | sed -r 's/.{10}$//g' | sed -r 's/.{4}//g')
        echo $Timepoint
        Condition=$(echo $RawData | rev | cut -d '/' -f1 | rev | sed -r 's/.{12}$//g' | sed 's/Fv//g')
        echo $Condition
        OutDir=qc_rna/RNAseq/fastqc/raw/$Condition/$Timepoint
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
        sbatch -p short $ProgDir/fastqc2.sh $RawData $OutDir
    done
```


```bash
    # Run fastq-mcf
    for RNADir in $(ls -d ../../../archives/2021_camb_general/20210128_Fvenenatum_CarbonRNAseq/rawdata); do
    FileNum=$(ls $RNADir/*_1.fq.gz | grep 'FvC4' | wc -l)
        for num in $(seq 1 $FileNum); do
            printf "\n"
            FileF=$(ls $RNADir/*1.fq.gz | grep 'FvC4' | head -n $num | tail -n1)
            FileR=$(ls $RNADir/*2.fq.gz | grep 'FvC4' | head -n $num | tail -n1)
            echo $FileF
            echo $FileR
            #Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1.fq.gz//g')
            #echo $Sample_Name
            Timepoint=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed -r 's/.{10}$//g' | sed -r 's/.{4}//g')
            echo $Timepoint
            Condition=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed -r 's/.{12}$//g' | sed 's/Fv//g')
            echo $Condition
            OutDir=qc_rna/RNAseq/Fvenenatum_CarbonRNAseq/$Condition/$Timepoint
            echo $OutDir
            IluminaAdapters=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc/illumina_full_adapters.fa
            ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
            sbatch -p himem $ProgDir/fastq-mcf_himem.sh $FileF $FileR $IluminaAdapters RNA $OutDir
        done
    done
```

```bash
# Run fastqc
    for RawData in $(ls qc_rna/RNAseq/Fvenenatum_CarbonRNAseq/*/*/*/*.fq.gz | grep 'FvC4'); do
        echo $RawData
        Timepoint=$(echo $RawData | rev | cut -d '/' -f3 | rev )
        echo $Timepoint
        Condition=$(echo $RawData | rev | cut -d '/' -f4 | rev )
        echo $Condition
        OutDir=qc_rna/RNAseq/fastqc/qc_rna/$Condition/$Timepoint
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
        sbatch -p short $ProgDir/fastqc2.sh $RawData $OutDir
    done
```

## Decontamination of rRNA reads in RNAseq data

```bash
    for RNADir in $(ls -d qc_rna/RNAseq/Fvenenatum_CarbonRNAseq/C*/T*); do
        FileNum=$(ls $RNADir/F/*_1_trim.fq.gz | wc -l)
        for num in $(seq 1 $FileNum); do
            printf "\n"
            FileF=$(ls $RNADir/F/*trim.fq.gz | head -n $num | tail -n1)
            FileR=$(ls $RNADir/R/*trim.fq.gz | head -n $num | tail -n1)
            echo $FileF
            echo $FileR
            Ref=/data/scratch/gomeza/prog/bbmap/ribokmers.fa.gz
            ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
            echo $RNADir
            Strain=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
            Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
            Condition=$(echo $RNADir | rev | cut -f2 -d '/' | rev)
            echo $Condition
            Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
            echo $Sample_Name
            echo $Timepoint
            echo $Strain
            sbatch -p himem $ProgDir/bbduk.sh $Ref "$RNADir"/cleaned/$Condition/$Timepoint/$Sample_Name $FileF $FileR $ProgDir $Strain
        done
    done
```

## Salmon 


```bash
conda activate salmon

for Transcriptome in $(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.cdna.fasta); do
Strain=$(echo $Transcriptome| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Transcriptome| rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for RNADir in $(ls -d ../../../data/scratch/gomeza/qc_rna/RNAseq/Fvenenatum_CarbonRNAseq/C*/T*/cleaned/C*/T*/* | grep 'FvC4'); do
FileNum=$(ls $RNADir/F/*_1_cleaned.fq.gz | wc -l)
for num in $(seq 1 $FileNum); do
printf "\n"
FileF=$(ls $RNADir/F/*cleaned.fq.gz | head -n $num | tail -n1)
FileR=$(ls $RNADir/R/*cleaned.fq.gz | head -n $num | tail -n1)
echo $FileF
echo $FileR
Prefix=$(echo $RNADir | rev | cut -f3 -d '/' | rev)
Timepoint=$(echo $RNADir | rev | cut -f2 -d '/' | rev)
Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_cleaned.fq.gz//g')
echo "$Prefix"
echo "$Timepoint"
echo "$Sample_Name"
OutDir=alignment/salmon/Fvenenatum_CarbonRNAseq/$Organism/$Strain/$Prefix/$Timepoint/$Sample_Name
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
sbatch -p himem $ProgDir/salmon.sh $Transcriptome $FileF $FileR $OutDir
done
done
done
```

