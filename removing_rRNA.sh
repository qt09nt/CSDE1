#this script is to remove rRNA contamination from the CSDE1 RIP sequencing
#using a software called RiboDetector
https://github.com/hzi-bifo/RiboDetector

#login to ARC
cd /work/yang_lab/1_Software

#activate conda environment
conda activate ribodetector

#try the CPU version of Ribodetector
## -e: Ensure which classification has high confidence for paired end reads.  norrna: output only high confident non-rRNAs, the rest are clasified as rRNAs;
ribodetector_cpu -t 20 \
-l 181 \
-i /work/yang_lab/queenie/Celf2-KI-Polysome-seq/fastq_files/Celf2-WT-mono-1_S2_R1_001.fastq \
-e norrna \ 
-o /work/yang_lab/queenie/Celf2-KI-Polysome-seq/ribodetector/Celf2-WI-mono-1.norrna.fastq

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

##############try GPU version of Ribodetector

#Request GPU core resource
#https://rcs.ucalgary.ca/How_to_request_an_interactive_GPU_on_ARC
salloc -N1 -n1 -c4 --mem=16GB --gres=gpu:1 -p gpu-v100 -t 1:00:00  

ribodetector -t 20 \
  -l 181 \
  -i /work/yang_lab/queenie/Celf2-KI-Polysome-seq/fastq_files/Celf2-WT-mono-1_S2_R1_001.fastq \
  -m 10 \
  -e norrna \ 
  --chunk_size 256 \
  -o /work/yang_lab/queenie/Celf2-KI-Polysome-seq/ribodetector/Celf2-WI-mono-1.norrna.fastq
