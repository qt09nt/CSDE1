https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit

#change directory to /work/yang_lab/1_Software/

#request compute node resources
salloc --mem=10G --time 01:00:00 --partition=single

#download SRA toolkit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.7/sratoolkit.3.0.7-centos_linux64.tar.gz

#extract toolkit
tar -xzf sratoolkit.3.0.7-centos_linux64.tar.gz

#get these samples
GSM1229978: M_P0_WT_0; Mus musculus; RNA-Seq (SRR980305)
SRX349397: GSM1229979: M_P0_WT_1; Mus musculus; RNA-Seq (SRR980306)
SRX349398: GSM1229980: M_P0_WT_2; Mus musculus; RNA-Seq (SRR980307)
