## Minion vs MiSeq. Feature comparison

All TFs identification and promoter analysis is done in the minion genome. This file contains commands used to compare the features annotated in the minion genome with the miseq genomes.

#### blast candidate TFs in miseq genome

```bash
mkdir Minion_vs_MiSeq
cd Minion_vs_MiSeq
```

```bash
# Extract g4106 g4107 and g6129 genes from fasta
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
$ProgDir/extract_from_fasta.py --fasta final_genes_appended_renamed.gene.fasta --headers head.txt > TFs.fasta
```

```bash
# Blast them on miseq genome
for Assembly in $(ls ../assembly/previous_versions/F.venenatum/WT/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f5 | rev)
echo "$Organism - $Strain"
Query=TFs.fasta
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
sbatch $ProgDir/blast_pipe.sh $Query dna $Assembly
done

# Blast them on miseq genes
for Assembly in $(ls ../gene_pred/codingquarry/F.venenatum/WT/final/final_genes_appended_renamed.gene.fasta); do
Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f5 | rev)
echo "$Organism - $Strain"
Query=TFs.fasta
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
sbatch $ProgDir/blast_pipe.sh $Query dna $Assembly
done

# output in dna/
```


```bash
# Extract g4106 g4107 and g6129 genes from pep file
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
$ProgDir/extract_from_fasta.py --fasta final_genes_appended_renamed.pep.fasta --headers head2.txt > TFs_v2.fasta
```

```bash
# Blast them on miseq genome
for Assembly in $(ls ../assembly/previous_versions/F.venenatum/WT/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f5 | rev)
echo "$Organism - $Strain"
Query=TFs_v2.fasta
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
sbatch $ProgDir/blast_pipe.sh $Query protein $Assembly
done
# Blast them on miseq genes
for Assembly in $(ls ../gene_pred/codingquarry/F.venenatum/WT/final/final_genes_appended_renamed.gene.fasta); do
Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f5 | rev)
echo "$Organism - $Strain"
Query=TFs_v2.fasta
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
sbatch $ProgDir/blast_pipe.sh $Query protein $Assembly
done

# output in prot/
```