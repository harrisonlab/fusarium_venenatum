# Co-expression analysis

### WGCNA

This analysis was done to compare the performance of WGCNA vs Cemitools


```bash
OutDir=GOMEZ/analysis/coexpression/WGCNA
#mkdir -p $OutDir
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Co-expression_analysis
gene_table=GOMEZ/analysis/coexpression/WGCNA/Carbon_data_vst_T_tab.txt
cat $gene_table | awk '{ print NF}'
column_start=2
column_end=146
/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/WGCNA_script.r --gene_table $gene_table --out_dir $OutDir --column_start $column_start --column_end $column_end
# Calculate softthreshold 
max_SFT=40
/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/softthreshold_power.r --out_dir $OutDir --max_SFT $max_SFT
# Cemitools scale-free topology model gave 12 score. WGCNA seems to be similar so I will use 12 too.
SFT=12
/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/create_network.r --out_dir $OutDir --sft $SFT --min_module_size 30 --merging_threshold 0.25
# Export network
/projects/software/R-3.6.1/bin/Rscript --vanilla $ProgDir/export_genes.r --out_dir $OutDir --unmerge Y
```


# Cemitools

Analysis performed on my local machine

```r
# This is the first run, so I tested multiple parameters and input data

# Input data and colData
reg<-read.table("vst_corrected.txt",header=TRUE,sep="\t")
reg <- reg %>% mutate_all(as.numeric)
raw<-read.table("FvenCarbon_RNAseq_design3.txt",header=T,sep="\t")
head(reg[,1:4])
cem <- cemitool(reg, raw)


TS1 <- as.data.frame.matrix("vst1_corrected.txt")


reg<-read.table("Carbon_data_vst_F.txt",header=TRUE)
reg<-read.table("vst1_corrected.txt",header=TRUE)
View(reg)
cem <- cemitool(reg)

vstT<-read.table("Carbon_data_vst_T.csv",header=TRUE,sep="\t")
cem <- cemitool(vstT, colData)

nmodules(cem)
generate_report(cem)
write_files(cem)
save_plots(cem, "all")

raw$Group <- paste0(raw$Class,'_', raw$Timepoint)


colData<-read.table("colData.txt",header=T,sep="\t")


vst_F does not work!!!

vst_T worked

I tried group and condition. both gave a thershold in 12 and not unique module in sucrose high

cem <- cemitool(vstT, colData,filter_pval=0.2) # it takes more time. One module more. No tri5
```