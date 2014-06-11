#!/bin/bash -eu

echo "*** Splitting BAM by chromosome ***"
echo "   " >> $LOGFILE
echo "*** Splitting BAM by chromosome ***" >> $LOGFILE

if [ $# -lt 3 ]
then
        echo "Usage: $0 <chr> <output> <bam>..."
        exit 1
fi

chr=$1
out=`cd \`dirname $2\`; pwd`/`basename $2`
shift 2

bams=''

for f in $*
do
	f=`cd \`dirname $f\`; pwd`/`basename $f`

	echo ">>> Extracting $chr from BAM: $f"
	echo ">>> Extracting $chr from BAM: $f" >> $LOGFILE
	o=${f/.bam/}.$chr.bam

	if [ $chr = 'UNK' -o $chr = 'chrU' -o $chr = 'U' ]
	then
		command="samtools view $f -f 12 -bo $o"
		samtools view $f -f 12 -bo $o
	else
		command="samtools view $f $chr -bo $o"
		samtools view $f $chr -bo $o
	fi
	echo ">>> $command"
	echo ">>> $command" >> $LOGFILE
	bams="$bams $o"
done

echo ">>> Merging $chr BAMs into $out"
echo ">>> Merging $chr BAMs into $out" >> $LOGFILE
if [ $# -gt 1 ]
then
	command="samtools merge $out $bams"
	samtools merge $out $bams
	echo ">>> $command"
	echo ">>> $command" >> $LOGFILE
	rm $bams
else
	mv $bams $out
fi

echo "*** Finished splitting BAM by chromosome ***"
