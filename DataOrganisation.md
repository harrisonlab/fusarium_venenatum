# Data migration to CropDiversity

This file contains all the data transfers done from our HPC to gruffalo. Only final data will be backed up

```bash
cd /data/scratch/gomeza
mkdir -p fusarium_venenatum/qc_rna/RNAseq/Carbon_dataset/
mv qc_rna/RNAseq/Fvenenatum_CarbonRNAseq/corrected/ fusarium_venenatum/qc_rna/RNAseq/Carbon_dataset/
scp -r fusarium_venenatum/qc_rna/RNAseq/Carbon_dataset/corrected agomez@gruffalo.cropdiversity.ac.uk:/home/agomez/projects/niab/agomez/fusarium_venenatum/4Backup/qc_rna/RNAseq/Carbon_dataset


```