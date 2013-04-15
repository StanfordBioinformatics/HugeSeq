#!/bin/bash -eu

echo "*** Generating alignment in BAM format ***"

if [ $# -lt 3 ]
then 
	echo "Usage: $0 <bam> <fastq1> <fastq2, '-' if none> [RG tag]"
	exit 1
fi

bam=`cd \`dirname $1\`; pwd`/`basename $1`
q1=`cd \`dirname $2\`; pwd`/`basename $2`
q2=`cd \`dirname $3\`; pwd`/`basename $3`

rgOpt=''
if [ $# -gt 3 ]
then
	rgOpt="-r $4"
fi

if [ $3 = '-' ]
then
	samse="bwa samse $rgOpt $REF $q1.sai $q1"
	echo ">> $samse"
	$samse | samtools view -bo $bam -S -
#	rm $q1.sai
else
	sampe="bwa sampe -P $rgOpt $REF $q1.sai $q2.sai $q1 $q2"
	echo ">> $sampe"
	$sampe | samtools view -bo $bam -S -
#	rm $q1.sai $q2.sai
fi

echo "*** Finished generating alignment in BAM format ***"
