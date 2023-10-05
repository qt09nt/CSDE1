#this script is for removing rRNA from the RNA seq files using Ribodetector

#https://github.com/hzi-bifo/RiboDetector

#Running it from a singularity container on ARC: The containerized version of ribodetector is available here: https://hub.docker.com/r/dawnmy/ribodetector
# instructions on converting a Docker container to the Singularity format here: https://rcs.ucalgary.ca/How_to_convert_a_Docker_container_to_an_Apptainer_container

#Converting the image to the Apptainer format
#Use apptainer build to pull and convert the Docker image
apptainer build ribodetector_image.sif docker://dawnmy/ribodetector:0.2.7

#move the Apptainer version of ribodetector to the software directory 
mv ribodetector_image.sif /work/yang_lab/1_Software/

#The command for launching ribodetector may require some modifications to enable GPU detection on the device. Additionally, 
#it's important to ensure that the container can access the yang_lab work directory. 

#Allocate resources and request GPU access:
[<user>@arc ~]$ salloc --mem=30G --time 01:00:00 -c 20 --partition=gpu-v100 --gres=gpu:1

#Change the directory to the yang_lab location
[<user>@fg1 ~]$ cd /work/yang_lab/1_Software/

#Run the application within the container using the apptainer command. Make sure to bind the ‘yang_lab’ directory to the container:

[<user>@fg1 1_Software]$ apptainer run --nv -B /work/yang_lab ribodetector_image.sif ribodetector -t 20 -l 181 \
          -i /work/yang_lab/queenie/Celf2-KI-Polysome-seq/fastq_files/Celf2-WT-mono-1_S2_R1_001.fastq -m 10 -e norrna \
          --chunk_size 256 -o /work/yang_lab/queenie/Celf2-KI-Polysome-seq/ribodetector/Celf2-WT-mono-1.norrna.fastq
