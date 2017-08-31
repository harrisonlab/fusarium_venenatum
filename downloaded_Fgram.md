# F. gramminearium

Informatics on the downloaded genome of Fusarium gramminearium

the assembly was downloaded from:
http://fungi.ensembl.org/Fusarium_graminearum/Info/Index?db=core
using assembly version RR1

Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/fusarium_venenatum

```bash
  cd /home/groups/harrisonlab/project_files/fusarium_venenatum
```

# Data download

Assemblies and gene models were downloaded for Fg from Emsembl at
http://fungi.ensembl.org/Fusarium_graminearum/Info/Index

these data were copied onto the cluster in the following locations:

```bash
  Organism=F.graminearum
  Strain=PH1
  mkdir -p assembly/external_group/$Organism/$Strain
  mkdir -p assembly/external_group/$Organism/$Strain/gff
```

Data was unzipped prior to further analysis:

```bash
  cd assembly/external_group/F.graminearum/PH1
  for Dir in $(ls -d *); do
    cd $Dir;
    gunzip *.gz;
    cd ../;
  done
  cd ../../../../
```

# Collecting assembly statistics


```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/external_group/F.graminearum/PH1/dna/Fusarium_graminearum.RR1.dna.toplevel.fa); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    echo "$Organism - $Strain"
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```


## Assessing the Gene space in the assembled genome:

```bash
for Assembly in $(ls assembly/external_group/F.graminearum/PH1/dna/Fusarium_graminearum.RR1.dna.toplevel.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
# BuscoDB="Fungal"
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```

```bash
	for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt); do  
		echo $File;
		cat $File | grep -e '(C)' -e 'Total';
	done
```
<!--
# identifying the % of the assembly that has been repeatmasked:

```bash

Unmasked=$(ls assembly/external_group/F.graminearum/PH1/dna/Fusarium_graminearum.RR1.dna.toplevel.fa)
cat $Unmasked | grep -v '>' | grep -o 'N' | wc -l
Hardmasked=$(ls assembly/external_group/F.graminearum/PH1/dna/Fusarium_graminearum.RR1.dna_rm.toplevel.fa)
cat $Hardmasked | grep -v '>' | grep -o 'N' | wc -l
``` -->


# Repeatmasking

Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The best assembly was used to perform repeatmasking

```bash
  for Assembly in $(ls assembly/external_group/F.graminearum/PH1/dna/Fusarium_graminearum.RR1.dna.toplevel.fa); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
    BestAss=assembly/spades/F.venenatum/WT/ncbi_edits/contigs_min_500bp_renamed.fasta
    OutDir=repeat_masked/F.venenatum/WT/illumina_assembly_ncbi
    qsub $ProgDir/rep_modeling.sh $BestAss $OutDir
    qsub $ProgDir/transposonPSI.sh $BestAss $OutDir
  done
```  

<!--
** % bases maked by repeatmasker: 4.75%**

** % bases masked by transposon psi: 4.19% **

The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and hardmasked files.

```bash

for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep -w 'WT' | grep 'ncbi'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa | grep -w 'WT' | grep 'ncbi'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```

```bash
for RepDir in $(ls -d repeat_masked/F.*/*/* | grep -w 'WT' | grep 'ncbi'); do
Strain=$(echo $RepDir | rev | cut -f2 -d '/' | rev)
Organism=$(echo $RepDir | rev | cut -f3 -d '/' | rev)  
RepMaskGff=$(ls $RepDir/*_contigs_hardmasked.gff)
TransPSIGff=$(ls $RepDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
printf "$Organism\t$Strain\n"
# printf "The number of bases masked by RepeatMasker:\t"
sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
# printf "The number of bases masked by TransposonPSI:\t"
sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
# printf "The total number of masked bases are:\t"
cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
echo
done
```

F.venenatum	WT
302604
144657
438768 -->


## Assessing the Gene space in predicted transcriptomes:

```bash
for Assembly in $(ls assembly/external_group/F.graminearum/PH1/cdna/Fusarium_graminearum.RR1.cdna.all.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
# BuscoDB="Fungal"
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/genes
qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```

```bash
	for File in $(ls gene_pred/busco/*/*/genes/*/short_summary_*.txt); do  
		echo $File;
		cat $File | grep -e '(C)' -e 'Total';
	done
```
