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

echo ">> Generating configuration file"
perl $BREAKDANCER/perl/bam2cfg.pl $bams > $p.cfg

echo ">> Performing read-pair mapping"
$BREAKDANCER/cpp/breakdancer_max $optO $p.cfg > $p.txt

echo ">> Converting output to GFF"
minsize=50
awkopt='
{
        size=$8; 
        start=($2<=$5?$2:$5);
        end=($5>=$2?$5:$2);
        tx="";
        if (size<0 && $7=="INS") 
                size=-size; 
        feat=$7; 
        if ($7=="DEL") 
                feat="Deletion"; 
        else if ($7=="INS") 
                feat="Insertion"; 
        else if ($7=="INV") 
                feat="Inversion";
        else if ($7=="ITX" || $7=="CTX") {
                feat="Translocation";
                start=$2;
                end=$2;
                size=-1;
                tx=";TX "$7";TCHR "$4";TSTART "$5
        }
        if (($7!="CTX" && $1==$4) || ($7=="CTX")) {
                if (size>='$minsize' || size==-1) 
                        print $1"\tBreakDancer\t"feat"\t"start"\t"end"\t.\t.\t.\tSize "size tx;
        }
}'

grep -v '\#' $p.txt | awk "$awkopt" | sort -k1,1 -k4n -k5n -u > $o

echo "*** Finished Calling SV using Read-Pair Mapping ***"
