#!/bin/bash -eu

echo "*** Aligning reads (finding SA coordiantes) ***"

if [ $# -lt 1 ]
then 
	echo "Usage: $0 <fastq> [num of threads]"
	exit 1
fi

q=30
f=`cd \`dirname $1\`; pwd`/`basename $1`

optI=""
if [[ $f =~ ".*_pf.f(ast)?q(\..*)?" ]]
then
	optI="-I"
fi

optT=""
if [ $# -gt 1 ]
then
	optT="-t $2"
fi

bwaaln="bwa aln $optI $optT -q $q $REF $f"

echo ">> $bwaaln"
$bwaaln > $f.sai

echo "*** Finished aligning reads (finding SA coordiantes) ***"
