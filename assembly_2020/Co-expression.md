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

