#this script is to remove rRNA contamination from the CSDE1 RIP sequencing
#using a software called RiboDetector
https://github.com/hzi-bifo/RiboDetector

#login to ARC
cd /work/yang_lab/1_Software

#get a CPU node resource
#salloc --mem=10G --time 05:00:00 --partition=single

#load conda modules
module load bioconda/conda3
module load biobuilds/2017.11

#activate conda environment
conda activate ribodetector

#try the CPU version of Ribodetector
## -e: Ensure which classification has high confidence for paired end reads.  norrna: output only high confident non-rRNAs, the rest are clasified as rRNAs;
ribodetector_cpu -t 20 \
-l 181 \
-i /work/yang_lab/queenie/Celf2-KI-Polysome-seq/fastq_files/Celf2-WT-mono-1_S2_R1_001.fastq \
-e norrna \ 
-o /work/yang_lab/queenie/Celf2-KI-Polysome-seq/ribodetector/Celf2-WT-mono-1.norrna.fastq

###error: 
#ValueError: This ORT build has ['AzureExecutionProvider', 'CPUExecutionProvider'] enabled. Since ORT 1.9, you are required to explicitly set the 
#providers parameter when instantiating InferenceSession. For example, onnxruntime.InferenceSession(..., providers=['AzureExecutionProvider', 'CPU  
#ExecutionProvider'], ...)


##### Try the GPU version of Ribodetector

#Install PyTorch from https://pytorch.org/get-started/locally/
pip3 install torch torchvision torchaudio

#allocate ARC gpu resource
#salloc --mem=10G --time 02:00:00 --partition=gpu-v100   ### this one produces ERROR:  error: Job submit/allocate failed: Job violates accounting/QOS policy (job submit limit, user's size and/or time limits)

salloc --mem=10G --time 01:00:00 --partition=bigmem

############## try GPU version of Ribodetector

#Request GPU core resource
#https://rcs.ucalgary.ca/How_to_request_an_interactive_GPU_on_ARC
salloc -N1 -n1 -c4 --mem=16GB --gres=gpu:1 -p gpu-v100 -t 1:00:00  

ribodetector -t 20 \
  -l 181 \
  -i /work/yang_lab/queenie/Celf2-KI-Polysome-seq/fastq_files/Celf2-WT-mono-1_S2_R1_001.fastq \
  -m 10 \
  -e norrna \ 
  --chunk_size 256 \
  -o /work/yang_lab/queenie/Celf2-KI-Polysome-seq/ribodetector/Celf2-WT-mono-1.norrna.fastq

### Running the GPU version of the command gives this error:
#RuntimeError: indices should be either on cpu or on the same device as the indexed tensor (cpu)


#find what version of cuda
conda list
#nvidia-cuda-runtime-cu11  11.7.99 

#check pytorch version installed
python -c "import torch; print(torch.__version__)"
#2.0.1+cu117

python -c 'import torch;print(torch.__version__);print(torch.version.cuda)'
#2.0.1+cu117
#11.7

#install mmcv 
pip install mmcv==2.0.0 -f https://download.openmmlab.com/mmcv/dist/cu117/torch2.0/index.html

#### re-running ribodetector GPU version after installing mmcv still gives this error: 
#RuntimeError: indices should be either on cpu or on the same device as the indexed tensor (cpu)


#check here for solutions:
https://github.com/open-mmlab/mmrotate/issues/511
https://github.com/AUTOMATIC1111/stable-diffusion-webui/issues/3958

#uninstall pytorch with 
conda uninstall pytorch
pip uninstall torch
pip uninstall torch

#try installing the same version of pytorch as what the Ribodetector software was tested on ie. from Ribodetector
#github it says "Our code was tested with pytorch v1.7, v1.7.1, v1.10.2."

#You should find the CUDA Version *highest CUDA version the installed driver supports on the top right corner of the comand's output. 
nvidia-smi

#output says CUDA Version: 11.4

#for installing other previous versions of Pytorch
#https://pytorch.org/get-started/previous-versions/
#v1.7.1
# CUDA 11.0
conda install pytorch==1.7.1 torchvision==0.8.2 torchaudio==0.7.2 cudatoolkit=11.0 -c pytorch
#error

#v1.7.0
# CUDA 11.0
conda install pytorch==1.7.0 torchvision==0.8.0 torchaudio==0.7.0 cudatoolkit=11.0 -c pytorch
##error

#try installing pytorch 1.10.1
conda install pytorch==1.10.1 torchvision==0.11.2 torchaudio==0.10.1 cudatoolkit=11.3 -c pytorch -c conda-forge
#ERROR:Solving environment: failed
# CondaValueError: Malformed version string '~': invalid character(s).

###try installing Pytorch v.2
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

#try running Ribodetector GPU version command again
ribodetector -t 20   -l 181   -i /work/yang_lab/queenie/Celf2-KI-Polysome-seq/fastq_files/Celf2-WT-mono-1_S2_R1_001.fastq   -m 10   -e norrna --chunk_size 256  -o /work/yang_lab/queenie/Celf2-KI-Polysome-seq/ribodetector/Celf2-WT-mono-1.norrna.fastq

#Traceback of TorchScript (most recent call last):
#  File "/home/queenie.tsang/.conda/envs/ribodetector/lib/python3.8/site-packages/ribodetector/model/model.py", l                                                      ine 118, in last_items
#    indices = sorted_last_indices(pack=pack)
#    if unsort and pack.unsorted_indices is not None:
#        indices = indices[pack.unsorted_indices]
#                  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ <--- HERE
#    return pack.data[indices]
#RuntimeError: indices should be either on cpu or on the same device as the indexed tensor (cpu)

#https://github.com/hzi-bifo/RiboDetector/issues/40
#https://github.com/hzi-bifo/RiboDetector/issues/34

#try installing CUDA if needed  - HPC support contacted for help
#https://developer.nvidia.com/cuda-downloads?target_os=Linux&target_arch=x86_64&Distribution=Rocky&target_version=8&target_type=runfile_local
#Linux x86_64 Rocky 8

#check linux distribution version
cat /etc/os-release

#check the version of CUDA with 
nvcc --version
##### try STAR

#https://github.com/alexdobin/STAR
# Get latest STAR source from releases
wget https://github.com/alexdobin/STAR/archive/2.7.11a.tar.gz
tar -xzf 2.7.11a.tar.gz
cd STAR-2.7.11a

export PATH=$PATH:/work/yang_lab/1_Software/STAR-2.7.11a/bin/Linux_x86_64/STAR:$PATH

#run STAR
./STAR

#error:
#./STAR: /lib64/libm.so.6: version `GLIBC_2.29' not found (required by ./STAR)
#./STAR: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.26' not found (required by ./STAR)

#check version of GLIBC 
ldd --version

#ldd (GNU libc) 2.28
#Copyright (C) 2018 Free Software Foundation, Inc.
#This is free software; see the source for copying conditions.  There is NO
#warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#Written by Roland McGrath and Ulrich Drepper.

#https://www.cyberithub.com/solved-glibc-2-29-not-found-error-on-ubuntu-linux/
#get glibc version 2.29
wget -c https://ftp.gnu.org/gnu/glibc/glibc-2.29.tar.gz

#extract glibc 
tar -zxvf glibc-2.29.tar.gz

cd glibc-2.29
mkdir glibc-build
cd glibc-build/

#configure the code for your local architecture by running below configure script.
../configure --prefix=/home/queenie.tsang/work/yang_lab/1_Software/STAR-2.7.11a/bin/Linux_x86_64/glibc-2.29/glibc-build

# compile the code using make command as shown below
make
# this "make" step takes several hours (around 4 hours)

#install using make install
make install

#https://github.com/alexdobin/STAR/issues/1484

##### glibc-2.29 is now installed but not sure how to link to this version instead of the default old version

########## Try running alignment with the older version of STAR_2.5.2b which is available from the module biobuilds/2017.11

#generate the genome index https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

#submit STAR build reference index as a batch job
cd /work/yang_lab/queenie
sbatch star.slurm


#when I run STAR with the transcriptome mouse.all.cds.fa file with the gtf file where ncRNA and rRNA genes have been removed with grep it gives this error:
##Fatal INPUT FILE error, no valid exon lines in the GTF file: mm10.filtered.out.rRNA.ncRNA.gtf
#Solution: check the formatting of the GTF file. Most likely cause is the difference in chromosome naming between GTF and FASTA file.

#apparently this error is because of : https://github.com/alexdobin/STAR/issues/421
#"You are using the transcriptome FASTA but genome GTF.
# do you want to map to the transcriptome (transcript sequences) or to the genome (chromosome sequences)? If the former, you do not need annotations GTF at all. If the latter, you would need to use genome FASTA (ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz) and genome GTF which you already have.
# Also, I would recommend using one of the latest versions of the GENCODE."

#so trying another way to remove rRNA where you construct a bowtie index using the rRNA database fasta file, then align to it, and exclude mapped reads (exclude reads
#that map to rRNA genes
#: https://www.biostars.org/p/377260/

#download the rRNA database file from SILVA https://www.arb-silva.de/download/arb-files/
#SILVA Release 138.1: Download the latest SILVA databases for ARB for small (16S/18S) and for large (23S/28S) subunit ribosomal RNAs.
#Released: 27.08.2020
#https://www.arb-silva.de/download/arb-files/
#	SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz	65 M
# download the large (23S/28S)ribosomal subunit database fasta file
wget 	https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz
gzip -d SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz

#download the SILVA aligned small (16S/18S, SSU) ribosomal subunit fasta file
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
gzip -d SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz

#Filtering rRNA from RNAseq data
#https://www.biostars.org/p/321714/#321727

#https://www.biostars.org/p/207311/
#### Try out SortmeRNA tool to remove rRNA
#https://github.com/sortmerna/sortmerna


### create new conda environment for sortmerna
#environment location: /home/queenie.tsang/.conda/envs/sortmerna_env

#activate sortmerna environment
conda activate sortmerna_env

#### 
conda install sortmerna
which sortmerna
#~/.conda/envs/sortmerna_env/bin/sortmerna

# test the installation
sortmerna --version
#SortMeRNA version 4.3.6

# view help
sortmerna -h

cd /work/yang_lab/queenie/csde1_RIP/

##### run sortmerna on interactive compute node
sortmerna --ref /work/yang_lab/queenie/SILVA_138.1_LSURef_NR99_tax_silva.fasta \
--ref /work/yang_lab/queenie/SILVA_138.1_SSURef_NR99_tax_silva.fasta \
--reads /work/yang_lab/queenie/csde1_RIP/Li35179_S1_R1_001.fastq 


################### Oct 2 2023 the bash script for sortmerna submitted Sept 29 ran  into an error 
#try running interactively

#takes a reallly long time to run; allocate more compute resources to task
salloc --partition=bigmem --time=16:0:0 --nodes=1 --ntasks=1 --cpus-per-task=40 --mem=0

module load bioconda/conda3
conda activate sortmerna_env

#use sortmerna Reference database 
#https://sortmerna.readthedocs.io/en/latest/databases.html
wget https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz
mkdir rRNA_databases_v4.3.6
tar -xvf database.tar.gz -C rRNA_databases_v4.3.6

#Purge kvdb directory prior to each new run:
rm -rf $HOME/sortmerna/run/kvdb/

# run sortmerna interactively with one of the reference databases from their github site
# the parameter –num-alignment INT = - Very fast for INT=1
#--fastx ->   Optional  Output aligned reads into FASTA/FASTQ file
sortmerna --ref /work/yang_lab/queenie/rRNA_databases_v4.3.6/smr_v4.3_default_db.fasta --reads /work/yang_lab/queenie/csde1_RIP/Li35179_S1_R1_001.fastq --fastx --aligned --other -num_alignments 1

#### it finished running but there is no output generated in /home/queenie.tsang/sortmerna/run/out; also there is an error: 
#[is_split_ready:726] found existing readfeed descriptor /home/queenie.tsang/sortmerna/run/readb/readfeed
#terminate called after throwing an instance of 'std::out_of_range'
#  what():  stoi
#Aborted (core dumped)

#commented this issue in the sortmerna github here:
#https://github.com/sortmerna/sortmerna/issues/379

#try a different file as input 
