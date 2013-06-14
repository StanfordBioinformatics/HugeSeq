#step1: SAMtools: Sort the bam file by name and output to sorted_by_name.bam
qsub -q extended -l h_vmem=6G -R y -b y /srv/gs1/projects/snyder/cuiping/app/samtools-0.1.14/samtools sort -n -m 2000000000 \
$inputBAM \
$outputBAM

#step2: PICARD: Convert sorted_by_name.bam to fastq 
qsub -q extended -l h_vmem=16G -R y -b y /usr/java/latest/bin/java -Xms13g -Xmx13g \
-jar /srv/gs1/projects/snyder/cuiping/app/picard-tools-1.32/SamToFastq.jar \
INPUT=/srv/gs1/projects/snyder/cuiping/data/schizophreniaTwins/SS6002854/sorted2.SS6002854.bam \
FASTQ=/srv/gs1/projects/snyder/cuiping/data/schizophreniaTwins/SS6002854/sorted2.SS6002854.read1.fastq \
SECOND_END_FASTQ=/srv/gs1/projects/snyder/cuiping/data/schizophreniaTwins/SS6002854/sorted2.SS6002854.read2.fastq \
NON_PF=true \
TMP_DIR=/srv/gs1/projects/snyder/cuiping/data/schizophreniaTwins/test/tmp



