#!/bin/bash -eu

echo "*** Performing VCF Merging ***"

if [ $# -lt 2 ]
then
	echo "Usage: $0 <output> <VCF file to concat>..."
	exit 1
fi

output=`cd \`dirname $1\`; pwd`/`basename $1`

shift

inputs=''
for i in $*
do
	i=`cd \`dirname $i\`; pwd`/`basename $i`
	inputs="$inputs $i"
done

echo ">> Concatenating VCFs into $output"
vcf-concat $inputs | vcf-sort > $output
