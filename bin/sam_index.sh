#!/bin/bash -eu

echo "*** Indexing BAM ***"

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
fi

echo "*** Finished indexing BAM ***"
