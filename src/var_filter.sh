#!/bin/bash -eu

echo "*** Performing Variant Filtration ***"

if [ $# -lt 1 ]
then 
	echo "Usage: $0 <vcf>"
	exit 1
fi

f=`cd \`dirname $1\`; pwd`/`basename $1`
o=${f/.vcf/}.filtered.vcf
filter="(AB ?: 0) > 0.75 || QUAL < 50.0 || DP > 360 || SB > -0.1 || MQ0>=4"

echo ">> Filtering variants using filter: $filter"
java -Xms5g -Xmx5g -Djava.io.tmpdir=$TMP -jar $GATK/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R $REF \
	-o $o \
	-B:variant,VCF $f \
	-et NO_ET \
	--clusterWindowSize 10 \
	--filterExpression "$filter" \
	--filterName "filter"

mv $o $f
rm -f $o.*

echo "*** Finished Variant Filtration ***"
