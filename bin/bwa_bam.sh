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

#bam=`cd \`dirname $1\`; pwd`/`basename $1`
#cores=$2
#PL=$3
#LB=$4
#SM=$5
#ID=$6

echo ">> BAM input"
echo ">> BAM input" >> $LOGFILE
command="samtools bam2fq $bam | bwa mem -CMp $optT $optRG $REF - | samtools view -Sbt $REF.fai -o ${bam/.bam/}.bwa.bam -"
samtools bam2fq $bam | bwa mem -CMp $optT $optRG $REF - | samtools view -Sbt $REF.fai -o ${bam/.bam/}.bwa.bam -
echo ">>> $command"
echo ">>> $command" >> $LOGFILE

echo "*** Finished aligning reads ***"


#if [[ ${f: -4} = ".bam" ]]
#then
#        echo ">> BAM input"
#        echo ">> BAM input" >> $LOGFILE
#        #samtools bam2fq $f | bwa mem -CMp $optT $optRG $REF - | samtools view -Sbt $REF.fai -o ${f/.bam/}.bwa.bam -
#        samtools bam2fq $f | bwa mem -CMp $optT $optRG $REF - | samtools view -Sbt $REF.fai -o ${f/.bam/}.bwa.bam -
#	
#elif [[ ${f: -6}==".fastq" || ${f: -9}==".fastq.gz" ]]
#then
#        q1=`cd \`dirname $1\`; pwd`/`basename $1`
#	
#	if [[ ${f: -6}==".fastq" ]]
#	then
#	exit
#		f=$(echo $f | sed -e "s/.fastq//g")
#		output="${f}bam"
#	fi
#	if [[ ${f: -9}==".fastq.gz" ]]
#	then
#	exit
#		f=$(echo $f | sed -e "s/.fastq.gz//g")
#	fi
#
#	command="bwa mem $REF $q1 $optT $optRG | samtools view -Sbt $REF.fai -o $f.bwa.bam -"
#	echo ">>> $command"
#	echo ">>> $command" >> $LOGFILE
#	$command
#        #echo $bwamem 
#fi
#
#echo "*** Finished aligning reads ***"
