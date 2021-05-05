OutDir=analysis/coexpression
#mkdir -p $OutDir
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Co-expression_analysis
gene_table=analysis/coexpression/TPM.txt
column_start=2
column_end=121
/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/WGCNA_script.r --gene_table $gene_table --out_dir $OutDir --column_start $column_start --column_end $column_end



OutDir=analysis/coexpression
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Co-expression_analysis
max_SFT=40
/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/softthreshold_power.r --out_dir $OutDir --max_SFT $max_SFT


OutDir=analysis/coexpression
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Co-expression_analysis
SFT=9
/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/create_network.r --out_dir $OutDir --sft $SFT --min_module_size 30 --merging_threshold 0.25


OutDir=analysis/coexpression_VAG/min50
mkdir -p $OutDir
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Co-expression_analysis
SFT=11
/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/create_network.r --out_dir $OutDir --sft $SFT --min_module_size 50 --merging_threshold 0.25


OutDir=analysis/coexpression
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Co-expression_analysis
/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/export_genes.r --out_dir $OutDir --unmerge Y

grep -E 'g6431.t1' analysis/coexpression/*/Genes_in_*
analysis/coexpression/merged_modules/Genes_in_yellow.txt:"g6431.t1"
analysis/coexpression/unmerged_modules/Genes_in_midnightblue.txt:"g6431.t1"

3851 genes in yellow
231 genes in midnightblue

OutDir=analysis/coexpression
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Co-expression_analysis
/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/export2cytoscape.r --out_dir $OutDir --module yellow

OutDir=analysis/coexpression
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Co-expression_analysis
/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/export2cytoscape.r --out_dir $OutDir --module midnightblue





OutDir=analysis/WGCNA
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Co-expression_analysis
gene_table=analysis/WGCNA/vst4WGCNA.txt
column_start=2
column_end=146
/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/WGCNA_script.r --gene_table $gene_table --out_dir $OutDir --column_start $column_start --column_end $column_end

OutDir=analysis/WGCNA
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Co-expression_analysis
max_SFT=40
/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/softthreshold_power.r --out_dir $OutDir --max_SFT $max_SFT


OutDir=analysis/WGCNA
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Co-expression_analysis
SFT=18
/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/create_network.r --out_dir $OutDir --sft $SFT --min_module_size 30 --merging_threshold 0.25


OutDir=analysis/WGCNA
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Co-expression_analysis
/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/export_genes.r --out_dir $OutDir --unmerge Y

grep -E 'g6431.t1' analysis/WGCNA/*/Genes_in_*
analysis/WGCNA/merged_modules/Genes_in_paleturquoise.txt:"g6431"
analysis/WGCNA/unmerged_modules/Genes_in_paleturquoise.txt:"g6431"



