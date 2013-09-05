#!/bin/bash -eu

echo "*** Splitting BAM by chromosome ***"

if [ $# -lt 3 ]
then
        echo "Usage: $0 <chr> <output> <bam>..."
        exit 1
fi

chr=$1
out=`cd \`dirname $2\`; pwd`/`basename $2`
echo $out
shift 2

bams=''

for f in $*
do
	f=`cd \`dirname $f\`; pwd`/`basename $f`

	echo ">>> Extracting $chr from BAM: $f"
	o=${f/.bam/}.$chr.bam

	if [ $chr = 'UNK' -o $chr = 'chrU' -o $chr = 'U' ]
	then
		samtools view $f -f 12 -bo $o
	else
		samtools view $f $chr -bo $o
	fi
	bams="$bams $o"
done

echo ">>> Merging $chr BAMs into $out"
if [ $# -gt 1 ]
then
	samtools merge $out $bams
	rm $bams
else
	mv $bams $out
fi

echo "*** Finished splitting BAM by chromosome ***"
