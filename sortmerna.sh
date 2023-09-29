#!/bin/bash
#<------------------------Request for Resources----------------------->
#SBATCH -J sortmerna
#SBATCH --partition=bigmem,cpu2017-bf05
#SBATCH --time=12:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80
#SBATCH --mem=0

module load biobuilds/2017.11

REF_SSU=/work/yang_lab/queenie/SILVA_138.1_SSURef_NR99_tax_silva.fasta
REF_LSU=/work/yang_lab/queenie/SILVA_138.1_LSURef_NR99_tax_silva.fasta
input=/work/yang_lab/csde1_RIP/*

# Run
i=0

for FILE in $input*fastq
do echo "sortmerna starting for"
echo $i
echo "Accessing this file"
echo $FILE

#clear the contents in the kvdb directory prior to each new run
rm -rf $HOME/sortmerna/run/kvdb/

##### run sortmerna 
sortmerna --ref $REF_SSU --ref $REF_LSU --reads /work/yang_lab/queenie/csde1_RIP/$FILE

echo "trial" $i "complete"
((i++))
echo "-------------"
done

