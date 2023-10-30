#!/bin/bash
#<------------------------Request for Resources----------------------->
#SBATCH -J ribodetector
#SBATCH --mem=30G
#SBATCH --time 01:00:00
#SBATCH -c 20
#SBATCH --partition=gpu-v100
#SBATCH --gres=gpu:1

#print second column (fastq file 2) name
cat samples.txt | while read line; do awk '{print $2}';  done

#save first column of line in samples.txt into the variable $file1
cat samples.txt | while read line; do file1=`awk '{print $1}'`; echo $file1;  done

cat samples.txt | while read line; do file1=`awk '{print $2}'`; echo $file1;  done
