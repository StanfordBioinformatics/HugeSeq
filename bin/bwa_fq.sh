#!/bin/bash -eu

echo "*** Aligning reads using BWA MEM algorithm ***"
echo "   " >> $LOGFILE
echo "*** Aligning reads using BWA MEM algorithm ***" > $LOGFILE

if [ $# -lt 1 ]
then 
	echo "Usage: $0 <fastq1> <fastq2> [num of threads] [RG tag]  or $0 <bam> [num of threads] [RG tag]"
	exit 1
fi

fq1=`cd \`dirname $1\`; pwd`/`basename $1`
fq2=`cd \`dirname $2\`; pwd`/`basename $2`

optRG=""
lastArg=${BASH_ARGV[0]}
if [[ $lastArg =~ "@RG" ]]
then
        optRG="-R $lastArg"
fi
optRG="-R @RG\tID:Default\tLB:Library\tPL:Illumina\tSM:SAMPLE"

optT=""
seclastArg=${@: -2:1}
optT="-t $seclastArg"

if [[ ${fq1: -6} == ".fastq" ]]
then
        fq=$(echo $fq1 | sed -e "s/.fastq//g")
        output="${fq}bam"
elif [[ ${fq1: -9} == ".fastq.gz" ]]
then
        fq=$(echo $fq1 | sed -e "s/.fastq.gz//g")
fi

command="bwa mem $REF $fq1 $fq2 $optT $optRG | samtools view -Sbt $REF.fai -o $fq.bwa.bam -"
echo ">>> $command"
echo ">>> $command" >> $LOGFILE
bwa mem $REF $fq1 $fq2 $optT $optRG | samtools view -Sbt $REF.fai -o $fq.bwa.bam -

echo "*** Finished aligning reads ***"

