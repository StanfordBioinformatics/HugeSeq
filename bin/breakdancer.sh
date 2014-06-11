#!/bin/bash -eu

echo "*** Calling SV using Read-Pair Mapping: $BREAKDANCER ***"

if [ $# -lt 2 ]
then
	echo "Usage: $0 <output> <bam>..."
	exit 1
fi

o=`cd \`dirname $1\`; pwd`/`basename $1`
p=${o/.gff/}
shift

bams=''
for i in $*
do
	bams="$bams `cd \`dirname $i\`; pwd`/`basename $i`"
done

optO=''
if [ $# -eq 1 ]
then
        if [[ "$1" =~ ".*chr[^.]*\..*" ]]
        then
                chr=`echo "$1" | sed 's/.*\(chr[^\.]*\)\..*/\1/'`
		# Segmentaion fault problem for -o, commented out temporarily
		#optO="-o $chr"
        fi
fi

echo ">> Generating configuration file $bams"
perl $BREAKDANCER/perl/bam2cfg.pl $bams > $p.cfg

echo ">> Performing read-pair mapping"
$BREAKDANCER/cpp/breakdancer_max $optO $p.cfg > $p.txt

echo ">> Converting output to GFF"
minsize=50
awkopt='{size=$8; if (size<0 && $7=="INS") size=-size; feat=$7; if ($7=="DEL") feat="Deletion"; else if ($7=="INS") feat="Insertion"; else if ($7=="INV") feat="Inversion"; if (feat!="ITX" && $1==$4 && size>='$minsize') print $1"\tBreakDancer\t"feat"\t"($2<=$5?$2:$5)"\t"($5>=$2?$5:$2)"\t"$9"\t.\t.\tSize "size"; nr.reads: "$10};'

#echo -e "#Chr\tProgram\tSV-type\t\tstart\tend\tscore\tstrand\tframe\tattributes" > $o
if [ -e "$o" ]
then
        rm $o
fi

grep -v '\#' $p.txt | awk "$awkopt" | sort -k1,1 -k4n -k5n -u >> $o

echo "*** Finished Calling SV using Read-Pair Mapping ***"
