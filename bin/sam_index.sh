#!/bin/bash -eu

echo "*** Indexing2 BAM ***"

if [ $# -lt 1 ]
then 
	echo "Usage: $0 <bam>"
	exit 1
fi

f=`cd \`dirname $1\`; pwd`/`basename $1`

if [ ! -e $f.bai ]
then
	echo ">>> BAM file $f is being indexed"
	samtools index $f
	echo ">>> Fixing BAI name"
	python $HUGESEQ_HOME/bin/fix_bai_name.py $f.bai
fi

echo "*** Finished indexing BAM ***"
