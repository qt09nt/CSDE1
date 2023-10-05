#!/bin/bash
#<------------------------Request for Resources----------------------->
#SBATCH -J ribodetector
#SBATCH --mem=30G
#SBATCH --time 01:00:00
#SBATCH -c 20
#SBATCH --partition=gpu-v100
#SBATCH --gres=gpu:1

cd /work/yang_lab/1_Software/

input=/work/yang_lab/queenie/csde1_RIP/fastq/*

# Run
i=0

for FILE in $input*fastq
do echo "ribodetector starting for"
echo $i
echo "Accessing this file"
echo $FILE

##### Run ribodetector GPU version
apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
          -i $FILE -m 12 -e norrna \
          --chunk_size 256 -o "$FILE.norrna.fastq"

#### move the output file to new folder
mv *.norrna.fastq /work/yang_lab/queenie/csde1_RIP/ribodetector/

"trial" $i "complete"
((i++))
echo "-------------"
done
