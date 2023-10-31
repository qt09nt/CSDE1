#!/bin/bash
#<------------------------Request for Resources----------------------->
#SBATCH -J ribodetector
#SBATCH --mem=30G
#SBATCH --time 01:00:00
#SBATCH -c 20
#SBATCH --partition=gpu-v100
#SBATCH --gres=gpu:1

#save first column of line in samples.txt into the variable $file1 save the second column (fastq 2 file) into file2 variable
#cat samples.txt | while read line; do file1=`awk '{print $1}'`; echo $file1;  done
#cat samples.txt | while read line; do file1=`awk '{print $2}'`; echo $file1;  done

cd /work/yang_lab/1_Software/

#input=/work/yang_lab/queenie/krausher_2014/fastq/
#output=/work/yang_lab/queenie/krausher_2014/ribodetector/

##### Run ribodetector GPU version

echo "T_P0_WT_0";

	apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
                -i "/work/yang_lab/queenie/krausher_2014/fastq/T_P0_WT_0_SRR980315_1.fastq" "/work/yang_lab/queenie/krausher_2014/fastq/T_P0_WT_0_SRR980315_2.fastq" \
                -m 12 \
                -e norrna \
                --chunk_size 256 \
                -o "/work/yang_lab/queenie/krausher_2014/ribodetector/T_P0_WT_0_SRR980315_1.norrna.fastq" "/work/yang_lab/queenie/krausher_2014/ribodetector/T_P0_WT_0_SRR980315_2.norrna.fastq"


echo "T_P0_WT_1";
  apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
                -i "/work/yang_lab/queenie/krausher_2014/fastq/T_P0_WT_1_SRR980316_1.fastq" "/work/yang_lab/queenie/krausher_2014/fastq/T_P0_WT_1_SRR980316_2.fastq" \
                -m 12 \
                -e norrna \
                --chunk_size 256 \
                -o "/work/yang_lab/queenie/krausher_2014/ribodetector/T_P0_WT_1_SRR980316_1.norrna.fastq" "/work/yang_lab/queenie/krausher_2014/ribodetector/T_P0_WT_1_SRR980316_2.norrna.fastq"

echo "T_P0_WT_2";
  apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
                -i "/work/yang_lab/queenie/krausher_2014/fastq/T_P0_WT_2_SRR980317_1.fastq" "/work/yang_lab/queenie/krausher_2014/fastq/T_P0_WT_2_SRR980317_2.fastq" \
                -m 12 \
                -e norrna \
                --chunk_size 256 \
                -o "/work/yang_lab/queenie/krausher_2014/ribodetector/T_P0_WT_2_SRR980317_1.norrna.fastq" "/work/yang_lab/queenie/krausher_2014/ribodetector/T_P0_WT_2_SRR980317_2.norrna.fastq"

echo "M_P0_WT_0";
  apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
                -i "/work/yang_lab/queenie/krausher_2014/fastq/M_P0_WT_0_SRR980305_1.fastq" "/work/yang_lab/queenie/krausher_2014/fastq/M_P0_WT_0_SRR980305_2.fastq" \
                -m 12 \
                -e norrna \
                --chunk_size 256 \
                -o "/work/yang_lab/queenie/krausher_2014/ribodetector/M_P0_WT_0_SRR980305_1.norrna.fastq" "/work/yang_lab/queenie/krausher_2014/ribodetector/M_P0_WT_0_SRR980305_2.norrna.fastq"

echo "M_P0_WT_1";
  apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
                -i "/work/yang_lab/queenie/krausher_2014/fastq/M_P0_WT_1_SRR980306_1.fastq" "/work/yang_lab/queenie/krausher_2014/fastq/M_P0_WT_1_SRR980306_2.fastq" \
                -m 12 \
                -e norrna \
                --chunk_size 256 \
                -o "/work/yang_lab/queenie/krausher_2014/ribodetector/M_P0_WT_1_SRR980306_1.norrna.fastq" "/work/yang_lab/queenie/krausher_2014/ribodetector/M_P0_WT_1_SRR980306_2.norrna.fastq"

echo "M_P0_WT_2";
  apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
                -i "/work/yang_lab/queenie/krausher_2014/fastq/M_P0_WT_2_SRR980307_1.fastq" "/work/yang_lab/queenie/krausher_2014/fastq/M_P0_WT_2_SRR980307_2.fastq" \
                -m 12 \
                -e norrna \
                --chunk_size 256 \
                -o "/work/yang_lab/queenie/krausher_2014/ribodetector/M_P0_WT_2_SRR980307_1.norrna.fastq" "/work/yang_lab/queenie/krausher_2014/ribodetector/M_P0_WT_2_SRR980307_2.norrna.fastq"

echo "P_P0_WT_0";
  apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
                -i "/work/yang_lab/queenie/krausher_2014/fastq/P_P0_WT_0_SRR980310_1.fastq" "/work/yang_lab/queenie/krausher_2014/fastq/P_P0_WT_0_SRR980310_2.fastq" \
                -m 12 \
                -e norrna \
                --chunk_size 256 \
                -o "/work/yang_lab/queenie/krausher_2014/ribodetector/P_P0_WT_0_SRR980310_1.norrna.fastq" "/work/yang_lab/queenie/krausher_2014/ribodetector/P_P0_WT_0_SRR980310_2.norrna.fastq"

echo "P_P0_WT_1";
  apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
                -i "/work/yang_lab/queenie/krausher_2014/fastq/P_P0_WT_1_SRR980311_1.fastq" "/work/yang_lab/queenie/krausher_2014/fastq/P_P0_WT_1_SRR980311_2.fastq" \
                -m 12 \
                -e norrna \
                --chunk_size 256 \
                -o "/work/yang_lab/queenie/krausher_2014/ribodetector/P_P0_WT_1_SRR980311_1.norrna.fastq" "/work/yang_lab/queenie/krausher_2014/ribodetector/P_P0_WT_1_SRR980311_2.norrna.fastq"

echo "P_P0_WT_2";
  apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
                -i "input=/work/yang_lab/queenie/krausher_2014/fastq/P_P0_WT_2_SRR980312_1.fastq" "input=/work/yang_lab/queenie/krausher_2014/fastq/P_P0_WT_2_SRR980312_2.fastq" \
                -m 12 \
                -e norrna \
                --chunk_size 256 \
                -o "/work/yang_lab/queenie/krausher_2014/ribodetector/P_P0_WT_2_SRR980312_1.norrna.fastq" "/work/yang_lab/queenie/krausher_2014/ribodetector/P_P0_WT_2_SRR980312_2.norrna.fastq"


