#!/bin/bash -eu

echo "*** Performing SNP Analysis using the GATK Unified Genotyper ***"

if [ $# -lt 2 ]
then 
	echo "Usage: $0 <output> <bam>..."
	exit 1
fi

o=`cd \`dirname $1\`; pwd`/`basename $1`
shift

f=''
for i in $*
do
        f="$f -I `cd \`dirname $i\`; pwd`/`basename $i`"
done

optL=''
if [ $# -eq 1 ]
then
        if [[ "$1" =~ ".*chr[^.]*\..*" ]]
        then
                chr=`echo "$1" | sed 's/.*\(chr[^\.]*\)\..*/\1/'`
                optL="-L $chr"
        fi
fi

echo ">>> Running the unified genotyper for SNP calling"
java -Xms5g -Xmx5g -Djava.io.tmpdir=$TMP -jar $GATK/GenomeAnalysisTK.jar \
	-T UnifiedGenotyper \
	-I $f \
	-R $REF \
	-D $SNP \
	-o $o $optL \
	-et NO_ET \
	-dcov 1000 \
	-A AlleleBalance \
	-A DepthOfCoverage \
	-A MappingQualityZero \
	-baq CALCULATE_AS_NECESSARY \
 	-stand_call_conf 30.0 \
	-stand_emit_conf 10.0

echo "*** Finished SNP Analysis using the GATK Unified Genotyper ***"
