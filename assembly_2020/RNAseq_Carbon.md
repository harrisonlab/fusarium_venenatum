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
# Run fastq-mcf
FileF=$(ls ../../../archives/2021_camb_general/20210128_Fvenenatum_CarbonRNAseq/rawdata/FvC0T0_1_1.fq.gz)
FileR=$(ls ../../../archives/2021_camb_general/20210128_Fvenenatum_CarbonRNAseq/rawdata/FvC0T0_1_2.fq.gz)
echo $FileF
echo $FileR
OutDir=qc_rna/RNAseq/Fvenenatum_CarbonRNAseq/FvC0/T0
IluminaAdapters=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc/illumina_full_adapters.fa
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
sbatch $ProgDir/fastq-mcf_himem.sh $FileF $FileR $IluminaAdapters RNA $OutDir
```

```bash
    # Run fastqc
    for RawData in $(ls qc_rna/48DD_experiment2020/WT53/T08/*/*.fq.gz | grep 'WT53_1\|WT53_2\|WT53_3\|WT53_4\|WT53_8') 
    do
        echo $RawData
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
        sbatch $ProgDir/fastqc.sh $RawData
    done
```

## Decontamination of rRNA reads i