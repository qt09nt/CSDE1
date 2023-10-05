#!/bin/bash
#<------------------------Request for Resources----------------------->
#SBATCH -J ribodetector
#SBATCH --mem=30G
#SBATCH --time 01:00:00 
#SBATCH -c 20 
#SBATCH --partition=gpu-v100 
#SBATCH --gres=gpu:1

cd /work/yang_lab/1_Software/

input=/work/yang_lab/queenie/csde1_RIP/*

# Run
i=0

for FILE in $input*fastq
do echo "sortmerna starting for"
echo $i
echo "Accessing this file"
echo $FILE

##### Run ribodetector GPU version
apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
          -i /work/yang_lab/queenie/Celf2-KI-Polysome-seq/fastq_files/${FILE}.fastq -m 12 -e norrna \
          --chunk_size 256 -o /work/yang_lab/queenie/Celf2-KI-Polysome-seq/ribodetector/${FILE}.norrna.fastq

"trial" $i "complete"
((i++))
echo "-------------"
done
