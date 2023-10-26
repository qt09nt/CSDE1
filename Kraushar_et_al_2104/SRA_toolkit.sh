https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit

#change directory to /work/yang_lab/1_Software/

#request compute node resources
salloc --mem=10G --time 01:00:00 --partition=single

#download SRA toolkit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.7/sratoolkit.3.0.7-centos_linux64.tar.gz

#extract toolkit
tar -xzf sratoolkit.3.0.7-centos_linux64.tar.gz

#configure SRA toolkit:
https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

#export the SRA bin directory onto the $PATH variable
export PATH=$PATH:/work/yang_lab/1_Software/sratoolkit.3.0.7-centos_linux64/bin

#get these samples
GSM1229978: M_P0_WT_0; Mus musculus; RNA-Seq (SRR980305)
SRX349397: GSM1229979: M_P0_WT_1; Mus musculus; RNA-Seq (SRR980306)
SRX349398: GSM1229980: M_P0_WT_2; Mus musculus; RNA-Seq (SRR980307)

SRX349401: GSM1229983: P_P0_WT_0; Mus musculus; RNA-Seq (SRR980310)
SRX349402: GSM1229984: P_P0_WT_1; Mus musculus; RNA-Seq (SRR980311)
SRX349403: GSM1229985: P_P0_WT_2; Mus musculus; RNA-Seq (SRR980312)

SRX349406: GSM1229988: T_P0_WT_0; Mus musculus; RNA-Seq (SRR980315)
SRX349407: GSM1229989: T_P0_WT_1; Mus musculus; RNA-Seq (SRR980316)
SRX349408: GSM1229990: T_P0_WT_2; Mus musculus; RNA-Seq (SRR980317)

#download the samples into this directory
/work/yang_lab/queenie/krausher_2014

#submit bash script for batch SRA download
sbatch auto_SRA.slurm

#rename fastq files with sample names

