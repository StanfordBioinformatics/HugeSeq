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

if [[ ${f: -4} = ".bam" ]]
then
        echo ">> BAM input"
        #samtools bam2fq $f | bwa mem -CMp $optT $optRG $REF - | samtools view -Sbt $REF.fai -o ${f/.bam/}.bwa.bam -
        samtools bam2fq $f | bwa mem -CMp $optT $optRG $REF - | samtools view -Sbt $REF.fai -o ${f/.bam/}.bwa.bam -
	
elif [[ ${f: -6}==".fastq" || ${f: -9}==".fastq.gz" ]]
then
        q1=`cd \`dirname $1\`; pwd`/`basename $1`
	
	if [[ ${f: -6}==".fastq" ]]
	then
	exit
		f=$(echo $f | sed -e "s/.fastq//g")
		output="${f}bam"
	fi
	if [[ ${f: -9}==".fastq.gz" ]]
	then
	exit
		f=$(echo $f | sed -e "s/.fastq.gz//g")
	fi

	bwa mem $REF $q1 $optT $optRG | samtools view -Sbt $REF.fai -o $f.bam -
        #echo $bwamem 
fi

echo "*** Finished aligning reads ***"
