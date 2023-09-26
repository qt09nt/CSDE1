#this script is to remove rRNA contamination from the CSDE1 RIP sequencing
#using a software called RiboDetector
https://github.com/hzi-bifo/RiboDetector

#login to ARC
cd /work/yang_lab/1_Software

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

