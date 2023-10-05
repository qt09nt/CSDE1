#!/bin/bash 

#SBATCH -J Celf2-KI-Polysome-seq
#SBATCH --mem=100G
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH -t 05:00:00
#SBATCH -c 16


#ASSIGN TOOLS
kallisto=/work/yang_lab/reza/tools/kallisto/kallisto
input=/work/yang_lab/queenie/csde1_RIP/ribodetector/*
output=/work/yang_lab/queenie/csde1_RIP/ribodetector/kallisto
reference=/work/yang_lab/queenie/mm10.refMrna


# Run
i=0

for FILE in $input*fastq 
do echo "trials starting for"
echo $i
echo "Accessing this file"
echo $FILE

$kallisto quant -i $reference -o $output$i".kallisto" --single -l 181 -s 20 --bias -b 100 -t 40 $FILE

echo "trial" $i "complete"
((i++))
echo "-------------"
done
done
#-----zip files---
echo "starting zip"
zip -r 20231005_csde1_ribodetector_kallisto_files.zip /work/yang_lab/queenie/csde1_RIP/
echo "zip done"
