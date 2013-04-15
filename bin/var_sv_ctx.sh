#!/bin/bash -eu

echo "*** Calling Transchromosomal SV using Read-Pair Mapping: $BREAKDANCER ***"

if [ $# -lt 2 ]
then
	echo "Usage: $0 <output> <bam>..."
	exit 1
fi

o=`cd \`dirname $1\`; pwd`/`basename $1`
p=${o/.gff/}
shift

bams=''
iOpt=''
for i in $*
do
	bam=`cd \`dirname $i\`; pwd`/`basename $i`
	bams="$bams $bam"
	iOpt="$iOpt I=$bam"
done

echo ">>> Merging input BAMs temporarily"
java -Xms5g -Xmx5g -jar $PICARD/MergeSamFiles.jar \
        $iOpt \
        O=${p}.tmp.bam \
        ASSUME_SORTED=true \
	VALIDATION_STRINGENCY=SILENT

echo ">> Generating configuration file"
perl $BREAKDANCER/perl/bam2cfg.pl ${p}.tmp.bam > $p.cfg

echo ">> Performing read-pair mapping"
$BREAKDANCER/cpp/breakdancer_max -t $p.cfg > $p.txt

echo ">> Converting output to GFF"
minsize=50
awkopt='
{
        if ($7=="CTX") {
                feat="Translocation";
                start=$2;
                end=$2;
                size=-1;
                tx=";TX "$7";TCHR "$4";TSTART "$5
                if (size>='$minsize' || size==-1) 
                        print $1"\tBreakDancer\t"feat"\t"start"\t"end"\t.\t.\t.\tSize "size tx;
        }
}'

grep -v '\#' $p.txt | awk "$awkopt" | sort -k1,1 -k4n -k5n -u > $o

rm -f ${p}.tmp.bam

echo "*** Finished Calling Transchromosomal SV using Read-Pair Mapping ***"
