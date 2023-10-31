###NOTE THIS SCRIPT STILL HAS ERRORS WHEN RUN 

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

# Run
i=0

##### Run ribodetector GPU version
        apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
                -i "/work/yang_lab/queenie/krausher_2014/fastq/T_P0_WT_0_SRR980315_1.fastq" "/work/yang_lab/queenie/krausher_2014/fastq/T_P0_WT_0_SRR980315_2.fastq" \
                -m 12 \
                -e norrna \
                --chunk_size 256 \
                -o "/work/yang_lab/queenie/krausher_2014/ribodetector/T_P0_WT_0_SRR980315_1.fastq.norrna.fastq" "/work/yang_lab/queenie/krausher_2014/ribodetector/T_P0_WT_0_SRR980315_2.fastq.norrna.fastq"



  apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
                -i "/work/yang_lab/queenie/krausher_2014/fastq/${file1}" "/work/yang_lab/queenie/krausher_2014/fastq/${file2}" \
                -m 12 \
                -e norrna \
                --chunk_size 256 \
                -o "/work/yang_lab/queenie/krausher_2014/ribodetector/${file1}.norrna.fastq" "/work/yang_lab/queenie/krausher_2014/ribodetector/{file2}.norrna.fastq"

  apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
                -i "/work/yang_lab/queenie/krausher_2014/fastq/${file1}" "/work/yang_lab/queenie/krausher_2014/fastq/${file2}" \
                -m 12 \
                -e norrna \
                --chunk_size 256 \
                -o "/work/yang_lab/queenie/krausher_2014/ribodetector/{file1}.norrna.fastq" "/work/yang_lab/queenie/krausher_2014/ribodetector/{file2}.norrna.fastq"

  apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
                -i "/work/yang_lab/queenie/krausher_2014/fastq/${file1}" "/work/yang_lab/queenie/krausher_2014/fastq/${file2}" \
                -m 12 \
                -e norrna \
                --chunk_size 256 \
                -o "/work/yang_lab/queenie/krausher_2014/ribodetector/{file1}.norrna.fastq" "/work/yang_lab/queenie/krausher_2014/ribodetector/{file2}.norrna.fastq"

  apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
                -i "/work/yang_lab/queenie/krausher_2014/fastq/${file1}" "/work/yang_lab/queenie/krausher_2014/fastq/${file2}" \
                -m 12 \
                -e norrna \
                --chunk_size 256 \
                -o "/work/yang_lab/queenie/krausher_2014/ribodetector/${file1}.norrna.fastq" "/work/yang_lab/queenie/krausher_2014/ribodetector/{file2}.norrna.fastq"

  apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
                -i "/work/yang_lab/queenie/krausher_2014/fastq/${file1}" "/work/yang_lab/queenie/krausher_2014/fastq/${file2}" \
                -m 12 \
                -e norrna \
                --chunk_size 256 \
                -o "/work/yang_lab/queenie/krausher_2014/ribodetector/{file1}.norrna.fastq" "/work/yang_lab/queenie/krausher_2014/ribodetector/{file2}.norrna.fastq"

  apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
                -i "/work/yang_lab/queenie/krausher_2014/fastq/${file1}" "/work/yang_lab/queenie/krausher_2014/fastq/${file2}" \
                -m 12 \
                -e norrna \
                --chunk_size 256 \
                -o "/work/yang_lab/queenie/krausher_2014/ribodetector/${file1}.norrna.fastq" "/work/yang_lab/queenie/krausher_2014/ribodetector/${file2}.norrna.fastq"

  apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
                -i "/work/yang_lab/queenie/krausher_2014/fastq/${file1}" "/work/yang_lab/queenie/krausher_2014/fastq/${file2}" \
                -m 12 \
                -e norrna \
                --chunk_size 256 \
                -o "/work/yang_lab/queenie/krausher_2014/ribodetector/${file1}.norrna.fastq" "/work/yang_lab/queenie/krausher_2014/ribodetector/${file2}.norrna.fastq"

  apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 100 \
                -i "${input}${file1}" "${input}${file2}" \
                -m 12 \
                -e norrna \
                --chunk_size 256 \
                -o "/work/yang_lab/queenie/krausher_2014/ribodetector/${file1}.norrna.fastq" "/work/yang_lab/queenie/krausher_2014/ribodetector/${file2}.norrna.fastq"

cd /work/yang_lab/1_Software/

echo "trial" $i "complete"
((i++))
echo "-------------"
done
