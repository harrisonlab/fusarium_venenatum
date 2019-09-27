#!/bin/bash
#SBATCH --partition=unlimited
#SBATCH --time=12:00:00
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=4


Assembly=$1
CassisTSV=$2
GeneID=$3
OutDir=$4

conda activate meme-v4

cassis \
  --annotation $CassisTSV \
  --genome $Assembly \
  --anchor $GeneID \
  --dir $OutDir/$GeneID \
  --mismatches 0 \
  -v \
  --prediction \
  --num-cpus 4 \
  | tee 2>&1 $OutDir/$GeneID/${GeneID}_log.txt
