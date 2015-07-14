#!/bin/bash -eu

echo "*** Calling CNV using Read-Depth Analysis: $CNVNATOR ***"

if [ $# -lt 2 ]
then
	echo "Usage: $0 <output> <root/bam...>"
	exit 1
fi

CNVNATOR=$CNVNATOR/cnvnator

binsize=100

o=`cd \`dirname $1\`; pwd`/`basename $1`
shift

bams=''
for i in $*
do
	bams="$bams `cd \`dirname $i\`; pwd`/`basename $i`"
done


if [ $# -eq 1 ]
then
        if [[ "$1" =~ ".*chr[^.]*\..*" ]]
        then
                chr=`echo "$1" | sed 's/.*\(chr[^\.]*\)\..*/\1/'`
		CNVNATOR="$CNVNATOR -chrom $chr"
        fi
fi

if [ -z "${1/*.root/}" ]
then
	echo ">> Processing root file: $1"
	p=${1/.root/}
else
	echo ">> Extracting read mapping from bam files: $CNVNATOR -root ${o/.gff/}.root -tree $bams"
	p=${o/.gff/}
	$CNVNATOR -root $p.root -tree $bams
fi

(
echo ">> Generating histogram"
$CNVNATOR -root $p.root -his $binsize -d `dirname $REF`
echo ">> Calculating statistics"
$CNVNATOR -root $p.root -stat $binsize
echo ">> RD signal partitioning"
$CNVNATOR -root $p.root -partition $binsize
echo ">> CNV calling"
$CNVNATOR -root $p.root -call $binsize | grep -v "\(WARN\)\|\(==\)" > $p.txt
) 2>&1 | grep -v bound | grep -v Zero | grep -v corrected

echo ">> Converting output to GFF"

minsize=50
awkopt='{feat=$4; if ($4=="deletion") feat="Deletion"; else if ($4=="duplication") feat="Duplication"; if ($5>='$minsize') print $1"\tCNVnator\t"feat"\t"$2"\t"$3"\t"$7" "$8" "$9" "$10" "$11"\t.\t.\tSize "$5";RD "$6};'

if [ -e "$o" ]
then
	rm $o
fi 
cut -f 1,2,3,4,5,6 $p.txt | sed 's/\(.*\)	\(.*\):\(.*\)-\(.*\)	\(.*\)	\(.*\)	\(.*\)	\(.*\)/\2	\3	\4	\1	\5	\6	p1: \7 | p2: \8/' | awk "$awkopt" | sort -k1,1 -k4n -k5n -u >> $o

echo "*** Finished CNV Calling using Read-Depth Analysis"
