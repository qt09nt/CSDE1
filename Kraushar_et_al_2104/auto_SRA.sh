#!/bin/bash

#<---------Resources------------>
#SBATCH -J sratoolkit
#SBATCH --mem=32G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -c 32

#add to path the SRA toolkit executable file 
export PATH=$PATH:$HOME/work/yang_lab/1_Software/sratoolkit.3.0.7-centos_linux64/bin

cd /work/yang_lab/queenie/krausher_2014

#retrieve the fastq-lite files using the accession numbers using prefetch command, then extract the data using fasterq-dump command
for sample in SRR980305 SRR980306 SRR980307 SRR980310 SRR980311 SRR980312 SRR980315 SRR980316 SRR980317;

do 
	
	prefetch $sample  && fasterq-dump $sample;
done;
