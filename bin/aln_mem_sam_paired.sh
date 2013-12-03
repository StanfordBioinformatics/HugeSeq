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
        optRG="-R $lastArg"
fi

optT=""
seclastArg=${@: -2:1}
optT="-t $seclastArg"

if [[ ${f: -6}==".fastq" || ${f: -9}==".fastq.gz" ]]
then
        q1=`cd \`dirname $1\`; pwd`/`basename $1`
	q2=`cd \`dirname $2\`; pwd`/`basename $2`

	if [[ ${f: -6} == ".fastq" ]]
	then
		f=$(echo $f | sed -e "s/.fastq//g")
		output="${f}bam"
	elif [[ ${f: -9} == ".fastq.gz" ]]
	then
		f=$(echo $f | sed -e "s/.fastq.gz//g")
	fi
	bwa mem $REF $q1 $q2 $optT $optRG | samtools view -Sbt $REF.fai -o $f.bwa.bam -
        #echo $bwamem 
fi

echo "*** Finished aligning reads ***"
