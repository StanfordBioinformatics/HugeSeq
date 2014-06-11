#!/bin/bash -eu

echo "*** Indexing2 BAM ***"
echo "   " >> $LOGFILE
echo "*** Indexing BAM ***" >> $LOGFILE

if [ $# -lt 1 ]
then 
	echo "Usage: $0 <bam>"
	exit 1
fi

f=`cd \`dirname $1\`; pwd`/`basename $1`

if [ ! -e $f.bai ]
then
	
	command="samtools index $f"
	echo ">>> BAM file $f is being indexed"
	echo ">>> BAM file $f is being indexed" >> $LOGFILE
	echo ">>> $command"
	$command
	
	echo ">>> Fixing BAI name"	
	echo ">>> Fixing BAI name" >> $LOGFILE	
	command="python $HUGESEQ_HOME/bin/fix_bai_name.py $f.bai"
	echo ">>> $command"
	echo ">>> $command" >> $LOGFILE
	$command
fi

echo "*** Finished indexing BAM ***"
