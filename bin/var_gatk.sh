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
java -Xmx9g -Xms9g -Djava.io.tmpdir=$TMP -jar $GATK/GenomeAnalysisTK.jar \
   -T UnifiedGenotyper \
   -R $REF \
   -I $f \
   --dbsnp $SNP \
   -o $o \
   -dcov 1000 \
   -stand_call_conf 30.0 \
   -stand_emit_conf 10.0 \
   -gt_mode DISCOVERY \
   --genotype_likelihoods_model BOTH \
   -nct 2

echo "*** Finished SNP Analysis using the GATK Unified Genotyper ***"
