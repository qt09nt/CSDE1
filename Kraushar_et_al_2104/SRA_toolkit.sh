https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit

#change directory to /work/yang_lab/1_Software/

#request compute node resources
salloc --mem=10G --time 01:00:00 --partition=single

#download SRA toolkit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.7/sratoolkit.3.0.7-centos_linux64.tar.gz

#extract toolkit
tar -xzf sratoolkit.3.0.7-centos_linux64.tar.gz
