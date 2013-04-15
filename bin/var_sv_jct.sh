#!/bin/bash -eu

echo "*** Calling SV using BreakSeq: $BREAKSEQ ***"

if [ $# -lt 2 ]
then
	echo "Usage: $0 <output> <bam>..."
	exit 1
fi

output=`cd \`dirname $1\`; pwd`/`basename $1`
shift

bams=''
for i in $*
do
        bams="$bams `cd \`dirname $i\`; pwd`/`basename $i`"
done


echo ">> Invoking the BreakSeq (Lite) program (Library: $BPLIB)"
$BREAKSEQ/bin/breakseq $output $bams

echo "*** Finished Calling SV using BreakSeq ***"
