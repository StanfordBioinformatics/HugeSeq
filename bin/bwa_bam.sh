#!/bin/bash -eu

echo "*** Aligning reads using BWA MEM algorithm ***"
echo "   " >> $LOGFILE
echo "*** Aligning reads using BWA MEM algorithm ***" >> $LOGFILE

if [ $# -lt 1 ]
then 
	echo "Usage: $0 <fastq1> <fastq2> [num of threads] [RG tag]  or $0 <bam> [num of threads] [RG tag]"
	exit 1
fi


bam=`cd \`dirname $1\`; pwd`/`basename $1`
optRG=""
lastArg=${BASH_ARGV[0]}
if [[ $lastArg =~ "@RG" ]]
then
        optRG="-R $lastArg"
fi

optT=""
seclastArg=${@: -2:1}
optT="-t $seclastArg"

echo ">> BAM input"
echo ">> BAM input" >> $LOGFILE
command="samtools bam2fq $bam | bwa mem -CMp $optT $optRG $REF - | samtools view -Sbt $REF.fai -o ${bam/.bam/}.bwa.bam -"
echo ">>> $command"
samtools bam2fq $bam | bwa mem -CMp $optT $optRG $REF - | samtools view -Sbt $REF.fai -o ${bam/.bam/}.bwa.bam -
echo ">>> $command" >> $LOGFILE
echo "*** Finished aligning reads ***"
