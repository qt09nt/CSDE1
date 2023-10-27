#MultiQC - Generate a single report which almagamates all fastqc results from all samples
#https://multiqc.info/docs/getting_started/quick_start

# Launch an interactive session:
$ salloc --mem=10G --time 01:00:00 --partition=single

Load bioconda

module load bioconda/conda3
cd /work/yang_lab/1_Software

Create new conda environment
conda create --name multiqc python=3.11

source activate multiqc

To install MultiQC within the multiqc virtual environment
pip install multiqc

#navigate to directory where the fastqc reports are located
cd /work/yang_lab/queenie/Celf2-KI-Polysome-seq/fastqc

#run multiqc
multiqc .

To download the multiqc reports (this is one single command) open a home terminal on local machine: 
rsync -v queenie.tsang@arc-dtn.ucalgary.ca:/work/yang_lab/queenie/Celf2-KI-Polysome-seq/fastqc/multiqc_report.html .

To download the multiqc data directory:
Rsync -axv queenie.tsang@arc-dtn.ucalgary.ca:/work/yang_lab/queenie/Celf2-KI-Polysome-seq/fastqc/multiqc_data .







