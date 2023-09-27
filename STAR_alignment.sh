#alignment of mouse genome reference file and gtf annotation file where rRNA and noncoding RNA removed
#submit as a batch job 



module load biobuilds/2017.11

#STAR build reference index
STAR --runMode genomeGenerate --genomeDir /work/yang_lab/queenie --genomeFastaFiles Mus_musculus.GRCm38.cds.all.fa --sjdbGTFfile mm10.filtered.out.rRNA.ncRNA.gtf
