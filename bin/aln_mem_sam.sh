#!/bin/bash -eu

echo "*** Aligning reads using BWA MEM algorithm ***"

if [ $# -lt 1 ]
then 
	echo "Usage: $0 <fastq1> <fastq2> [num of threads] [RG tag]  or $0 <bam> [num of threads] [RG tag]"
	exit 1
fi


f=`cd \`dirname $1\`; pwd`/`basename $1`


optRG=""
lastArg=${BASH_ARGV[0]}
if [[ $lastArg =~ "@RG" ]]
then
        optRG="-r $lastArg"
fi
echo $optRG

optT=""
seclastArg=${@: -2:1}
if [[ $seclastArg =~ "^[0-9]+$" ]]
then
        optT="-t $seclastArg"
fi
echo $seclastArg


if [[ ${f: -4} = ".bam" ]]
then
        echo ">> BAM input"
        bwamem="`samtools bam2fq $f | bwa mem -CMp $optT $optRG $REF - | samtools view -Sbt $REF.fai -o ${f/.bam/}.bwa.bam -`"
        echo $bwamem
elif [[ ${f: -6}==".fastq" || ${f: -9}==".fastq.gz" ]]
then
        q1=`cd \`dirname $2\`; pwd`/`basename $2`
	q2=`cd \`dirname $3\`; pwd`/`basename $3`
        echo ">> FASTQ input"
	bwamem="\`bwa mem $REF $q1 $q2 | samtools view -Sbt $REF.fai -o ${f/.bam/}.bwa.bam -\`"
        echo ">> FASTQ input"
        echo $bwamem 
fi


echo "*** Finished aligning reads ***"
