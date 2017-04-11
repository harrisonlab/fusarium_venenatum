# F. gramminearium

Informatics on the downloaded genome of Fusarium gramminearium

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
  CurDir=$PWD
  cd assembly/external_group/$Organism/$Strain/dna
  gunzip Fusarium_graminearum.RR.dna.toplevel.fa.gz
  cd ../gff
  gunzip Fusarium_graminearum.RR.34.gff3.gz
  cd ../pep
  gunzip Fusarium_graminearum.RR.pep.all.fa.gz
  cd ../cds
  gunzip Fusarium_graminearum.RR.cds.all.fa.gz
  cd $CurDir
```

# Collecting assembly statistics


```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/external_group/F.graminearum/PH1/dna/Fusarium_graminearum.RR.dna.toplevel.fa); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    echo "$Organism - $Strain"
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```


## Assessing the Gene space in the assembled genome:

```bash
	for Assembly in $(ls assembly/external_group/F.graminearum/PH1/dna/Fusarium_graminearum.RR.dna.toplevel.fa); do
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
	for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt); do  
		echo $File;
		cat $File | grep -e '(C)' -e 'Total';
	done
```


## Assessing the Gene space in predicted transcriptomes:

```bash
	for Assembly in $(ls assembly/external_group/F.graminearum/PH1/cdna/Fusarium_graminearum.RR.cdna.all.fa); do
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
