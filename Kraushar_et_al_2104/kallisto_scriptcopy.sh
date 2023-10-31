#!/bin/bash

#SBATCH -J krausher_kallisto
#SBATCH --mem=100G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 05:00:00
#SBATCH -c 16


#ASSIGN TOOLS
kallisto=/work/yang_lab/reza/tools/kallisto/kallisto
input=/work/yang_lab/queenie/krausher_2014/ribodetector/
output=/work/yang_lab/queenie/krausher_2014/kallisto/
reference=/work/yang_lab/queenie/mm10.refMrna


# Run

#-s, --sd=DOUBLE               Estimated standard deviation of fragment length
#                              (default: -l, -s values are estimated from paired
#                               end data, but are required when using --single)
#-t, --threads=INT             Number of threads to use (default: 1)
#-b, --bootstrap-samples=INT   Number of bootstrap samples (default: 0)

echo "starting kallisto for M_P0_WT_ \n"
$kallisto quant -i $reference -o $output$i".kallisto" -s 20 -l 101 --bias -b 100 ${input}M_P0_WT_0_SRR980305_1.norrna.fastq ${input}M_P0_WT_0_SRR980305_2.norrna.fastq -t 40

echo "starting trial for M_P0_WT_1 \n"
$kallisto quant -i $reference -o $output$i".kallisto" -s 20 -l 101 --bias -b 100 ${input}M_P0_WT_1_SRR980306_1.norrna.fastq ${input}M_P0_WT_1_SRR980306_2.norrna.fastq -t 40

echo "starting trial for M_P0_WT_2 \n"
$kallisto quant -i $reference -o $output$i".kallisto" -s 20 -l 101 --bias -b 100 ${input}M_P0_WT_2_SRR980307_1.norrna.fastq ${input}M_P0_WT_2_SRR980307_2.norrna.fastq -t 40

echo "starting trial for P_P0_WT_0 \n"
$kallisto quant -i $reference -o $output$i".kallisto" -s 20 -l 101 --bias -b 100 ${input}P_P0_WT_0_SRR980310_1.norrna.fastq ${input}P_P0_WT_0_SRR980310_2.norrna.fastq -t 40

echo "starting trial for T_P0_WT_0 \n"
$kallisto quant -i $reference -o $output$i".kallisto" -s 20 -l 101 --bias -b 100 ${input}T_P0_WT_0_SRR980315_1.norrna.fastq ${input}T_P0_WT_0_SRR980315_2.norrna.fastq -t 40

echo "starting trial for T_P0_WT_1 \n"
$kallisto quant -i $reference -o $output$i".kallisto" -s 20 -l 101 --bias -b 100 ${input}T_P0_WT_1_SRR980316_1.norrna.fastq ${input}T_P0_WT_1_SRR980316_2.norrna.fastq -t 40

echo "starting trial forT_P0_WT_2 \n"
$kallisto quant -i $reference -o $output$i".kallisto" -s 20 -l 101 --bias -b 100 ${input}T_P0_WT_2_SRR980317_1.norrna.fastq ${input}T_P0_WT_2_SRR980317_2.norrna.fastq -t 40




~
~
~
~
~
~
~
~
~
"kallisto_scriptcopy.sh" 49L, 2244C                                                                                                                                                                                  45,67         All
