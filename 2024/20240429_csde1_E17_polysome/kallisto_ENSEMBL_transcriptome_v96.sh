#!/bin/bash

#SBATCH -J CSDE1_polysome_kallisto_ENSEMBL_v96
#SBATCH --mem=100G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 05:00:00
#SBATCH --cpus-per-task=40

#ASSIGN TOOLS
kallisto=/work/yang_lab/reza/tools/kallisto/kallisto
input=/work/yang_lab/2_Raw_data_yanglab/240429_csde1_E17_polysome/Data/fastq
output=/work/yang_lab/queenie/csde1_E17_polysome_20240429
reference=/work/yang_lab/queenie/mus_musculus/transcriptome.idx

$kallisto quant -i $reference -o $output"Csde1-WT1-total_S1.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Li47790-Csde1-WT1-total_S1_R1_001.fastq.gz
$kallisto quant -i $reference -o $output"Csde1-WT1-mono_S2.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Li47791-Csde1-WT1-mono_S2_R1_001.fastq.gz
$kallisto quant -i $reference -o $output"Csde1-WT1-poly_S3.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Li47792-Csde1-WT1-poly_S3_R1_001.fastq.gz
$kallisto quant -i $reference -o $output"Csde1-Het1-total_S4.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Li47793-Csde1-Het1-total_S4_R1_001.fastq.gz
$kallisto quant -i $reference -o $output"Csde1-Het1-mono_S5.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Li47794-Csde1-Het1-mono_S5_R1_001.fastq.gz
$kallisto quant -i $reference -o $output"Csde1-Het1-poly_S6.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Li47795-Csde1-Het1-poly_S6_R1_001.fastq.gz
$kallisto quant -i $reference -o $output"Csde1-WT3-total_S7.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Li47796-Csde1-WT3-total_S7_R1_001.fastq.gz
$kallisto quant -i $reference -o $output"Csde1-WT3-mono_S8.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Li47797-Csde1-WT3-mono_S8_R1_001.fastq.gz
$kallisto quant -i $reference -o $output"Csde1-WT3-poly_S9.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Li47798-Csde1-WT3-poly_S9_R1_001.fastq.gz
$kallisto quant -i $reference -o $output"Csde1-Het3-total_S10.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Li47799-Csde1-Het3-total_S10_R1_001.fastq.gz
$kallisto quant -i $reference -o $output"Csde1-Het3-mono_S11.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Li47800-Csde1-Het3-mono_S11_R1_001.fastq.gz
$kallisto quant -i $reference -o $output"Csde1-Het3-poly_S12.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Li47801-Csde1-Het3-poly_S12_R1_001.fastq.gz
$kallisto quant -i $reference -o $output"Csde1-WT4-total_S13.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Li47802-Csde1-WT4-total_S13_R1_001.fastq.gz
$kallisto quant -i $reference -o $output"Csde1-WT4-mono_S14.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Li47803-Csde1-WT4-mono_S14_R1_001.fastq.gz
$kallisto quant -i $reference -o $output"Csde1-WT4-poly_S15.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Li47804-Csde1-WT4-poly_S15_R1_001.fastq.gz
$kallisto quant -i $reference -o $output"Csde1-Het4-total_S16.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Li47805-Csde1-Het4-total_S16_R1_001.fastq.gz
$kallisto quant -i $reference -o $output"Csde1-Het4-mono_S17.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Li47806-Csde1-Het4-mono_S17_R1_001.fastq.gz
$kallisto quant -i $reference -o $output"Csde1-Het4-poly_S18.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Li47807-Csde1-Het4-poly_S18_R1_001.fastq.gz
$kallisto quant -i $reference -o $output"CHGI-NTC_S19.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Li47808-CHGI-NTC_S19_R1_001.fastq.gz
$kallisto quant -i $reference -o $output"Undetermined_S0.kallisto" --single -l 76 -s 20 --bias -b 100 -t 40 $input/Undetermined_S0_R1_001.fastq.gz
