#!/bin/bash -eu

echo "*** Performing VCF Merging ***"

if [ $# -lt 2 ]
then
	echo "Usage: $0 <output> <VCF file to merge>..."
	exit 1
fi

output=`cd \`dirname $1\`; pwd`/`basename $1`

shift

inputs=''
for i in $*
do
	i=`cd \`dirname $i\`; pwd`/`basename $i`
	echo ">> Indexing VCF: $i"
	bgzip -f $i
	tabix -f -p vcf $i.gz
	inputs="$inputs $i.gz"
done

echo ">> Merging VCFs into $output"
vcf-merge $inputs > $output
