# Codon Usage

Commands used to perform codon usage analysis of Fv




maria's documentation to prepare cds files:
https://github.com/harrisonlab/popgen/blob/master/codon/run_codonw.sh

maria's commands to run codonw
https://github.com/harrisonlab/popgen/blob/master/codon/codonw.sh


## Analysis of all genes:

```bash
cd /home/groups/harrisonlab/project_files/fusarium_venenatum
CDS=$(ls $PWD/gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.cds.fasta)
mkdir analysis/codon_usage
cd analysis/codon_usage
cp -s $CDS .
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/codon
$ProgDir/codonw.sh $(basename $CDS)
cd /home/groups/harrisonlab/project_files/fusarium_venenatum
ls /home/groups/harrisonlab/project_files/fusarium_venenatum/analysis/codon_usage/final_genes_appended_renamed.cds/codon.coa
```

Make a codon usage table to be usable by optimizer:
http://genomes.urv.es/OPTIMIZER/tutorial.php

```bash
CodonFile=$(ls analysis/codon_usage/final_genes_appended_renamed.cds/final_genes_appended_renamed.cds_pass.cutot)
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_venenatum/codon_usage
$ProgDir/codonw2optimizer.py --codonw $CodonFile

```

## Analysis of genes with high expression:


<!-- ```bash
cd /home/groups/harrisonlab/project_files/fusarium_venenatum
CDS=$(ls $PWD/gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.cds.fasta)
OutDir=analysis/codon_usage/fpkm_5
mkdir $OutDir
Fpkm=$(ls /home/deakig/projects/quorn/DGE/Quorn_fpkm.txt)
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_venenatum/codon_usage
$ProgDir/filter_fpkm.py --fpkm $Fpkm --threshold 5 > $OutDir/quorn_fpkm_5.tsv
cat $OutDir/quorn_fpkm_5.tsv | cut -f1 > $OutDir/quorn_fpkm_5_gene_IDs.txt
# A single transcript is taken from each gene
cat $OutDir/quorn_fpkm_5_gene_IDs.txt | tail -n+2 | sed "s/$/.t1/g" > $OutDir/quorn_fpkm_5_IDs.txt
CDS=$(ls $PWD/gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.cds.fasta)
# cat $CDS | grep -w -f $OutDir/quorn_fpkm_5_IDs.txt > $OutDir/quorn_fpkm_5_IDs.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $CDS --headers $OutDir/quorn_fpkm_5_IDs.txt > $OutDir/quorn_fpkm_5.fa
``` -->

```bash
OutDir=analysis/codon_usage/fpkm_148
mkdir $OutDir
Fpkm=$(ls /home/deakig/projects/quorn/DGE/Quorn_fpkm.txt)
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_venenatum/codon_usage
$ProgDir/filter_fpkm.py --fpkm $Fpkm --min 5 > $OutDir/quorn_fpkm_148.tsv
cat $OutDir/quorn_fpkm_148.tsv | cut -f1 > $OutDir/quorn_fpkm_148_gene_IDs.txt
# A single transcript is taken from each gene
cat $OutDir/quorn_fpkm_148_gene_IDs.txt | tail -n+2 | sed "s/$/.t1/g" > $OutDir/quorn_fpkm_148_IDs.txt
CDS=$(ls $PWD/gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.cds.fasta)
# cat $CDS | grep -w -f $OutDir/quorn_fpkm_148_IDs.txt > $OutDir/quorn_fpkm_148_IDs.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $CDS --headers $OutDir/quorn_fpkm_148_IDs.txt > $OutDir/quorn_fpkm_148.fa
```

```bash
AnnotTab=$(ls gene_pred/annotation/F.venenatum/WT/WT_annotation_ncbi.tsv)
cat $AnnotTab | grep -e 'IPR000477' -e 'IPR004875' -e 'IPR025476' -e 'IPR012337' -e 'transpos' > $OutDir/putative_transposon_IDs.txt

cat $AnnotTab | cut -f1,17,18 | grep -P -v "\t\w" | cut -f1 > $OutDir/no_annotation_IDs.txt

cat $OutDir/quorn_fpkm_148_IDs.txt | grep -v -f $OutDir/putative_transposon_IDs.txt | grep -v -f $OutDir/no_annotation_IDs.txt > $OutDir/quorn_fpkm_148_IDs_filtered.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $CDS --headers $OutDir/quorn_fpkm_148_IDs_filtered.txt > $OutDir/quorn_fpkm_148_filtered.fa
```


```bash
cd /home/groups/harrisonlab/project_files/fusarium_venenatum
OutDir=analysis/codon_usage/fpkm_148
cd $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/codon
$ProgDir/codonw.sh quorn_fpkm_148_filtered.fa
cd /home/groups/harrisonlab/project_files/fusarium_venenatum
ls $OutDir/*/quorn_fpkm_148_filtered_pass.cutot
```

Remove genes with transposon annotations:
IPR000477 - reverse transriptase
IPR004875 - DDE superfamily endonuclease
IPR025476 - Helitron helicase-like domain
IPR012337 - Ribonuclease H-like superfamily

Remove genes with no annotations:




Make a codon usage table to be usable by optimizer:
http://genomes.urv.es/OPTIMIZER/tutorial.php

```bash
CodonFile=$(ls $OutDir/*/quorn_fpkm_148_filtered_pass.cutot)
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_venenatum/codon_usage
$ProgDir/codonw2optimizer.py --codonw $CodonFile
```

```
UUU: 0.76; UCU: 1.28; UAU: 0.70; UGU: 0.78; UUC: 1.24; UCC: 1.08; UAC: 1.30; UGC: 1.22; UUA: 0.26; UCA: 0.97; UAA: 1.36; UGA: 0.80; UUG: 0.92; UCG: 0.80; UAG: 0.85; UGG: 1.00; CUU: 1.45; CCU: 1.48; CAU: 0.87; CGU: 1.16; CUC: 1.77; CCC: 1.19; CAC: 1.13; CGC: 1.44; CUA: 0.52; CCA: 0.86; CAA: 0.88; CGA: 1.59; CUG: 1.08; CCG: 0.47; CAG: 1.12; CGG: 0.47; AUU: 1.13; ACU: 1.09; AAU: 0.60; AGU: 0.70; AUC: 1.59; ACC: 1.27; AAC: 1.40; AGC: 1.17; AUA: 0.29; ACA: 1.07; AAA: 0.47; AGA: 0.78; AUG: 1.00; ACG: 0.56; AAG: 1.53; AGG: 0.55; GUU: 1.31; GCU: 1.40; GAU: 0.98; GGU: 1.36; GUC: 1.59; GCC: 1.34; GAC: 1.02; GGC: 1.43; GUA: 0.39; GCA: 0.74; GAA: 0.75; GGA: 0.92; GUG: 0.71; GCG: 0.52; GAG: 1.25; GGG: 0.30
```



## Genes with low expression

```bash
cd /home/groups/harrisonlab/project_files/fusarium_venenatum
OutDir=analysis/codon_usage/fpkm_1-10
mkdir $OutDir
Fpkm=$(ls /home/deakig/projects/quorn/DGE/Quorn_fpkm.txt)
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_venenatum/codon_usage
$ProgDir/filter_fpkm.py --fpkm $Fpkm --min 1 --max 10 > $OutDir/quorn_fpkm_1-10.tsv
cat $OutDir/quorn_fpkm_1-10.tsv | cut -f1 > $OutDir/quorn_fpkm_1-10_gene_IDs.txt
# A single transcript is taken from each gene
cat $OutDir/quorn_fpkm_1-10_gene_IDs.txt | tail -n+2 | sed "s/$/.t1/g" > $OutDir/quorn_fpkm_1-10_IDs.txt
# CDS=$(ls $PWD/gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.cds.fasta)
# # cat $CDS | grep -w -f $OutDir/quorn_fpkm_1-10_IDs.txt > $OutDir/quorn_fpkm_1-10_IDs.txt
# ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
# $ProgDir/extract_from_fasta.py --fasta $CDS --headers $OutDir/quorn_fpkm_1-10_IDs.txt > $OutDir/quorn_fpkm_1-10.fa
```

Remove genes with transposon annotations:
IPR000477 - reverse transriptase
IPR004875 - DDE superfamily endonuclease
IPR025476 - Helitron helicase-like domain
IPR012337 - Ribonuclease H-like superfamily

Remove genes with no annotations:

```bash
cd /home/groups/harrisonlab/project_files/fusarium_venenatum
AnnotTab=$(ls gene_pred/annotation/F.venenatum/WT/WT_annotation_ncbi.tsv)
cat $AnnotTab | grep -e 'IPR000477' -e 'IPR004875' -e 'IPR025476' -e 'IPR012337' -e 'transpos' > $OutDir/putative_transposon_IDs.txt

cat $AnnotTab | cut -f1,17,18 | grep -P -v "\t\w" | cut -f1 > $OutDir/no_annotation_IDs.txt

cat $OutDir/quorn_fpkm_1-10_IDs.txt | grep -v -f $OutDir/putative_transposon_IDs.txt | grep -v -f $OutDir/no_annotation_IDs.txt > $OutDir/quorn_fpkm_1-10_IDs_filtered.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $CDS --headers $OutDir/quorn_fpkm_1-10_IDs_filtered.txt > $OutDir/quorn_fpkm_1-10_filtered.fa
```


```bash
cd /home/groups/harrisonlab/project_files/fusarium_venenatum
OutDir=analysis/codon_usage/fpkm_1-10
cd $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/codon
$ProgDir/codonw.sh quorn_fpkm_1-10_filtered.fa
cd /home/groups/harrisonlab/project_files/fusarium_venenatum
ls $OutDir/*/quorn_fpkm_1-10_filtered_pass.cutot
```



Make a codon usage table to be usable by optimizer:
http://genomes.urv.es/OPTIMIZER/tutorial.php

```bash
CodonFile=$(ls $OutDir/*/quorn_fpkm_1-10_filtered_pass.cutot)
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_venenatum/codon_usage
$ProgDir/codonw2optimizer.py --codonw $CodonFile
```

```
UUU: 0.90; UCU: 1.22; UAU: 0.87; UGU: 0.89; UUC: 1.10; UCC: 0.89; UAC: 1.13; UGC: 1.11; UUA: 0.41; UCA: 1.15; UAA: 0.92; UGA: 1.23; UUG: 1.04; UCG: 0.87; UAG: 0.85; UGG: 1.00; CUU: 1.41; CCU: 1.34; CAU: 1.03; CGU: 1.03; CUC: 1.46; CCC: 0.95; CAC: 0.97; CGC: 1.21; CUA: 0.64; CCA: 1.12; CAA: 1.03; CGA: 1.44; CUG: 1.04; CCG: 0.59; CAG: 0.97; CGG: 0.60; AUU: 1.09; ACU: 1.03; AAU: 0.79; AGU: 0.80; AUC: 1.40; ACC: 1.01; AAC: 1.21; AGC: 1.07; AUA: 0.50; ACA: 1.25; AAA: 0.72; AGA: 1.02; AUG: 1.00; ACG: 0.70; AAG: 1.28; AGG: 0.69; GUU: 1.26; GCU: 1.28; GAU: 1.04; GGU: 1.17; GUC: 1.38; GCC: 1.10; GAC: 0.96; GGC: 1.31; GUA: 0.55; GCA: 1.02; GAA: 0.94; GGA: 1.06; GUG: 0.81; GCG: 0.59; GAG: 1.06; GGG: 0.47
```


## Optimise codons for given proteins

### GFP

```bash
OutDir=analysis/codon_usage/by_expression/optimised
mkdir -p $OutDir
printf ">GFP\nMVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSKLSKDPNEKRDHMVLLEFVTAAGITLGMDELYK\n" > $OutDir/GFP_unoptimised.aa

printf ">human_optimised_GFP\nATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCAAGCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG\n" > $OutDir/GFP_unoptimised.cds
printf ">Fv_GFP_optimizer_guided_random_fw\nATGGTGAGTAAAGGTGAAGAACTTTTCACAGGAGTGGTGCCAATCTTGGTAGAATTGGACGGAGACGTCAATGGACACAAGTTTTCAGTTTCGGGCGAAGGCGAGGGCGACGCAACGTACGGTAAGCTCACTCTTAAATTTATTTGCACAACGGGTAAACTCCCCGTGCCTTGGCCTACACTCGTAACTACGTTGACGTACGGTGTCCAATGCTTCTCTCGCTACCCCGATCATATGAAGCAGCACGATTTTTTTAAGTCTGCCATGCCTGAAGGCTACGTCCAAGAACGCACAATTTTCTTTAAGGACGATGGTAACTACAAGACACGAGCCGAGGTTAAGTTCGAAGGTGACACACTAGTGAACCGCATAGAACTCAAGGGCATTGATTTTAAGGAAGACGGAAATATCCTCGGGCATAAGCTTGAGTATAATTACAATTCTCATAACGTGTATATCATGGCTGATAAGCAAAAGAACGGCATCAAGGTCAACTTTAAGATTCGACACAATATCGAAGATGGCTCTGTACAGCTTGCTGACCATTATCAGCAAAACACACCAATCGGTGACGGTCCTGTGCTTCTACCTGATAACCACTACTTATCCACTCAAAGTAAACTGTCGAAGGACCCCAACGAAAAACGGGACCACATGGTCCTGCTGGAGTTTGTGACTGCAGCGGGCATCACACTTGGCATGGATGAGCTCTACAAG\n" >> $OutDir/GFP_unoptimised.cds

printf 'UUU: 0.76; UCU: 1.28; UAU: 0.70; UGU: 0.78; UUC: 1.24; UCC: 1.08; UAC: 1.30; UGC: 1.22; UUA: 0.26; UCA: 0.97; UAA: 1.36; UGA: 0.80; UUG: 0.92; UCG: 0.80; UAG: 0.85; UGG: 1.00; CUU: 1.45; CCU: 1.48; CAU: 0.87; CGU: 1.16; CUC: 1.77; CCC: 1.19; CAC: 1.13; CGC: 1.44; CUA: 0.52; CCA: 0.86; CAA: 0.88; CGA: 1.59; CUG: 1.08; CCG: 0.47; CAG: 1.12; CGG: 0.47; AUU: 1.13; ACU: 1.09; AAU: 0.60; AGU: 0.70; AUC: 1.59; ACC: 1.27; AAC: 1.40; AGC: 1.17; AUA: 0.29; ACA: 1.07; AAA: 0.47; AGA: 0.78; AUG: 1.00; ACG: 0.56; AAG: 1.53; AGG: 0.55; GUU: 1.31; GCU: 1.40; GAU: 0.98; GGU: 1.36; GUC: 1.59; GCC: 1.34; GAC: 1.02; GGC: 1.43; GUA: 0.39; GCA: 0.74; GAA: 0.75; GGA: 0.92; GUG: 0.71; GCG: 0.52; GAG: 1.25; GGG: 0.30' > $OutDir/Fv_codon_weights.cds

ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_clocks/codon_optimisation
$ProgDir/score_codons.py --fasta_aa $OutDir/GFP_unoptimised.aa  --fasta_cds $OutDir/GFP_unoptimised.cds --codon_table $OutDir/Fv_codon_weights.cds --prefix $OutDir/GFP
```

```
randomised codons:
AUGGUAAGUAAGGGAGAAGAACUCUUCACAGGGGUGGUGCCCAUAUUAGUUGAGUUAGACGGUGAUGUAAACGGACACAAGUUCUCCGUAUCGGGGGAGGGGGAGGGGGAUGCUACUUACGGCAAGCUAACUCUGAAGUUCAUAUGUACCACCGGGAAGCUCCCUGUACCGUGGCCCACAUUGGUAACCACACUAACAUAUGGUGUACAAUGCUUCUCGCGCUACCCGGAUCAUAUGAAGCAGCACGAUUUUUUCAAAUCCGCCAUGCCUGAGGGGUACGUGCAGGAGCGGACAAUUUUCUUCAAGGACGACGGCAACUACAAAACUCGCGCUGAGGUAAAGUUCGAGGGAGAUACAUUGGUAAACAGAAUAGAACUUAAGGGGAUAGAUUUCAAGGAGGAUGGGAACAUACUCGGACACAAGCUAGAGUACAAUUACAACAGCCACAACGUUUACAUUAUGGCUGACAAGCAGAAGAACGGGAUCAAGGUCAAUUUCAAGAUACGUCACAACAUAGAAGAUGGAUCGGUGCAACUUGCUGAUCAUUAUCAGCAGAACACCCCAAUAGGCGAUGGGCCGGUGCUGCUGCCCGACAACCAUUACCUCUCGACGCAGUCGAAGCUAUCUAAGGACCCAAACGAGAAGCGUGACCAUAUGGUACUUCUUGAGUUCGUGACGGCGGCAGGGAUAACACUCGGCAUGGACGAGCUUUACAAG
0.81 of 1.0
optimum sequence:
AUGGUCUCUAAGGGCGAGGAGCUCUUCACCGGCGUCGUCCCUAUCCUCGUCGAGCUCGACGGCGACGUCAACGGCCACAAGUUCUCUGUCUCUGGCGAGGGCGAGGGCGACGCUACCUACGGCAAGCUCACCCUCAAGUUCAUCUGCACCACCGGCAAGCUCCCUGUCCCUUGGCCUACCCUCGUCACCACCCUCACCUACGGCGUCCAGUGCUUCUCUCGAUACCCUGACCACAUGAAGCAGCACGACUUCUUCAAGUCUGCUAUGCCUGAGGGCUACGUCCAGGAGCGAACCAUCUUCUUCAAGGACGACGGCAACUACAAGACCCGAGCUGAGGUCAAGUUCGAGGGCGACACCCUCGUCAACCGAAUCGAGCUCAAGGGCAUCGACUUCAAGGAGGACGGCAACAUCCUCGGCCACAAGCUCGAGUACAACUACAACUCUCACAACGUCUACAUCAUGGCUGACAAGCAGAAGAACGGCAUCAAGGUCAACUUCAAGAUCCGACACAACAUCGAGGACGGCUCUGUCCAGCUCGCUGACCACUACCAGCAGAACACCCCUAUCGGCGACGGCCCUGUCCUCCUCCCUGACAACCACUACCUCUCUACCCAGUCUAAGCUCUCUAAGGACCCUAACGAGAAGCGAGACCACAUGGUCCUCCUCGAGUUCGUCACCGCUGCUGGCAUCACCCUCGGCAUGGACGAGCUCUACAAG
1.0 of 1.0
worst sequence:
AUGGUAAGUAAAGGGGAAGAAUUAUUUACGGGGGUAGUACCGAUAUUAGUAGAAUUAGAUGGGGAUGUAAAUGGGCAUAAAUUUAGUGUAAGUGGGGAAGGGGAAGGGGAUGCGACGUAUGGGAAAUUAACGUUAAAAUUUAUAUGUACGACGGGGAAAUUACCGGUACCGUGGCCGACGUUAGUAACGACGUUAACGUAUGGGGUACAAUGUUUUAGUCGGUAUCCGGAUCAUAUGAAACAACAUGAUUUUUUUAAAAGUGCGAUGCCGGAAGGGUAUGUACAAGAACGGACGAUAUUUUUUAAAGAUGAUGGGAAUUAUAAAACGCGGGCGGAAGUAAAAUUUGAAGGGGAUACGUUAGUAAAUCGGAUAGAAUUAAAAGGGAUAGAUUUUAAAGAAGAUGGGAAUAUAUUAGGGCAUAAAUUAGAAUAUAAUUAUAAUAGUCAUAAUGUAUAUAUAAUGGCGGAUAAACAAAAAAAUGGGAUAAAAGUAAAUUUUAAAAUACGGCAUAAUAUAGAAGAUGGGAGUGUACAAUUAGCGGAUCAUUAUCAACAAAAUACGCCGAUAGGGGAUGGGCCGGUAUUAUUACCGGAUAAUCAUUAUUUAAGUACGCAAAGUAAAUUAAGUAAAGAUCCGAAUGAAAAACGGGAUCAUAUGGUAUUAUUAGAAUUUGUAACGGCGGCGGGGAUAACGUUAGGGAUGGAUGAAUUAUAUAAA
0.49 of 1.0
midpoint sequence:
0.75
AUGGUGUCCAAGGGGGAGGAGUUGUUCACGGGAGUGGUGCCAAUACUCGUUGAGCUAGACGGCGACGUAAAUGGCCAUAAAUUUUCGGUGAGUGGGGAAGGGGAGGGAGACGCCACAUACGGAAAGCUAACACUGAAGUUCAUAUGCACCACCGGAAAGCUUCCAGUGCCGUGGCCGACUCUGGUCACGACACUCACAUACGGGGUACAAUGCUUUUCACGGUAUCCCGAUCAUAUGAAGCAGCAUGACUUCUUUAAAUCAGCAAUGCCGGAGGGCUACGUACAGGAACGAACGAUCUUCUUUAAGGAUGACGGGAACUACAAAACACGGGCGGAGGUAAAGUUUGAGGGGGAUACGCUGGUAAACAGGAUUGAACUAAAAGGCAUCGACUUCAAAGAAGACGGCAACAUCCUCGGGCAUAAGCUUGAAUAUAACUAUAACUCUCACAACGUGUACAUAAUGGCCGAUAAACAAAAGAACGGAAUAAAGGUAAACUUCAAAAUACGACAUAAUAUAGAAGACGGAAGCGUGCAGUUGGCCGACCAUUAUCAACAAAAUACGCCCAUAGGGGACGGGCCAGUGCUACUCCCUGAUAACCAUUACCUUUCGACGCAAUCAAAGCUUUCGAAGGACCCAAACGAGAAACGGGAUCACAUGGUCCUACUAGAAUUUGUGACGGCCGCAGGGAUCACCCUCGGGAUGGAUGAGCUAUACAAG
Score of the pre-optimised sequence:
human_optimised_GFP
0.93 of 1.0
Fv_GFP_optimizer_guided_random_fw
0.82 of 1.0
```

open an R session:

```bash
R
```

```R
GFP_1000_scores <- read.table("~/cluster_mount/groups/harrisonlab/project_files/fusarium_venenatum/analysis/codon_usage/by_expression/optimised/GFP_1000_scores.tsv", quote="\"", comment.char="")
library("ggplot2")
p <- ggplot(data=GFP_1000_scores, aes(GFP_1000_scores[,1])) + geom_histogram(bins = '50')
p <- p + xlim(0.49, 1.0)
p <- p + xlab("Total codon score") + ylab("Count")
p <- p + xlab("Total codon score") + ylab("Count")

p <- p + geom_vline(xintercept = 0.82, na.rm = FALSE, show.legend = NA)
p <- p + geom_vline(xintercept = 0.93, na.rm = FALSE, show.legend = NA)
p
```
