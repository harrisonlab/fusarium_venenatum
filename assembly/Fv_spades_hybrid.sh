#!/bin/bash

#Assemble contigs using SPAdes

#$ -S /bin/bash
#$ -cwd
#$ -pe smp 24
#$ -l virtual_free=15.7G
#$ -l h=blacklace11.blacklace


Usage="subSPAdes_2lib_pacbio.sh <Pacbio.fq> <F1_read.fa> <R1_read.fa> <F2_read.fa> <R2_read.fa> <output_directory> [<coverage_cutoff>]"
echo "$Usage"

M1=$1
M2=$2
F1=$3
R1=$4
F2=$5
R2=$6
F3=$7
R3=$8
OutDir=$9
Cutoff=20

CurPath=$PWD
WorkDir="$TMPDIR"

MinION_Read1=$(basename $M1)
MinION_Read2=$(basename $M2)
F1_Read=$(basename $F1)
R1_Read=$(basename $R1)
F2_Read=$(basename $F2)
R2_Read=$(basename $R2)
F3_Read=$(basename $F3)
R3_Read=$(basename $R3)

cp $CurPath/$M1 $WorkDir/$MinION_Read1
cp $CurPath/$M2 $WorkDir/$MinION_Read2
cp $CurPath/$F1 $WorkDir/$F1_Read
cp $CurPath/$R1 $WorkDir/$R1_Read
cp $CurPath/$F2 $WorkDir/$F2_Read
cp $CurPath/$R2 $WorkDir/$R2_Read
cp $CurPath/$F3 $WorkDir/$F3_Read
cp $CurPath/$R3 $WorkDir/$R3_Read

echo  "Running SPADES with the following inputs:"
echo "MinION reads = $M1\t$M2"
echo "F1_Read = $F1"
echo "R1_Read = $R1"
echo "F2_Read = $F2"
echo "R2_Read = $R2"
echo "F3_Read = $F3"
echo "R3_Read = $R3"
echo "Output directory will be: $CurPath/$OutDir"
echo "Coverage cutoff set to $Cutoff"


spades.py \
    -k 21,33,55,77,99,127 \
    -m 375 \
    --phred-offset 33 \
    --careful \
    --nanopore $WorkDir/$MinION_Read1 \
    --nanopore $WorkDir/$MinION_Read2 \
    --pe1-1 $WorkDir/$F1_Read \
    --pe1-2 $WorkDir/$R1_Read \
    --pe2-1 $WorkDir/$F2_Read \
    --pe2-2 $WorkDir/$R2_Read \
    --pe3-1 $WorkDir/$F3_Read \
    --pe3-2 $WorkDir/$R3_Read \
    -t 24  \
    -o $WorkDir/. \
    --cov-cutoff "$Cutoff"

rm $WorkDir/$MinION_Read1
rm $WorkDir/$MinION_Read2
rm $WorkDir/$F1_Read
rm $WorkDir/$R1_Read
rm $WorkDir/$F2_Read
rm $WorkDir/$R2_Read
rm $WorkDir/$F3_Read
rm $WorkDir/$R3_Read
mkdir -p $CurPath/$OutDir
cp -r $WorkDir/* $CurPath/$OutDir/.
echo "files copied"
