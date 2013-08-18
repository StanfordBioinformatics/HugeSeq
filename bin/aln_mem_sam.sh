#!/bin/bash -eu

echo "*** Aligning reads using BWA MEM algorithm ***"

if [ $# -lt 1 ]
then 
	echo "Usage: $0 <fastq1> <fastq2> [num of threads] [RG tag]  or $0 <bam> [num of threads] [RG tag]"
	exit 1
fi


f=`cd \`dirname $1\`; pwd`/`basename $1`
echo $f
optRG=""
lastArg=${BASH_ARGV[0]}
if [[ $lastArg =~ "@RG" ]]
then
        optRG="-R $lastArg"
fi
echo $optRG

optT=""
seclastArg=${@: -2:1}
#if [[ $seclastArg =~ "^[0-9]+$" ]]
#then
        optT="-t $seclastArg"
#fi
#echo $seclastArg
echo $optT

if [[ ${f: -4} = ".bam" ]]
then
        echo ">> BAM input"
        bwamem="`samtools bam2fq $f | bwa mem -CMp $optT $optRG $REF - | samtools view -Sbt $REF.fai -o ${f/.bam/}.bwa.bam -`"
        echo $bwamem
elif [[ ${f: -6}==".fastq" || ${f: -9}==".fastq.gz" ]]
then
        q1=`cd \`dirname $1\`; pwd`/`basename $1`
	#q2=`cd \`dirname $2\`; pwd`/`basename $2`
	echo $q1
	
	if [[ ${f: -6}==".fastq" ]]
	then
		f=$(echo $f | sed -e "s/.fastq//g")
		output="${f}bwa.bam"
	fi
	if [[ ${f: -9}==".fastq.gz" ]]
	then
		f=$(echo $f | sed -e "s/.fastq.gz//g")
	fi
	echo $f
	echo $optRG

	exit()	
	#bwamem="`bwa mem $REF $q1 $q2 $optT $optRG | samtools view -Sbt $REF.fai -o ${f/.bam/}.bwa.bam -`"
	bwamem="`bwa mem $REF $q1 $optT $optRG | samtools view -Sbt $REF.fai -o ${f/.bam/}.bwa.bam -`"
        echo $bwamem 
fi


echo "*** Finished aligning reads ***"
