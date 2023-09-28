#alignment of mouse genome reference file and gtf annotation file where rRNA and noncoding RNA removed
#submit as a batch job 

#!/bin/bash
#<------------------------Request for Resources----------------------->
#SBATCH --partition=bigmem, cpu2017-bf05
#SBATCH --time=6:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -cpus-per-task=8
#SBATCH --mem=0

module load biobuilds/2017.11

REF="Mus_musculus.GRCm38.cds.all.fa"
ANNOTATION="mm10.filtered.out.rRNA.ncRNA.gtf"

#STAR build reference index
STAR --runMode genomeGenerate --genomeDir /work/yang_lab/queenie --genomeFastaFiles $REF --sjdbGTFfile $ANNOTATION
