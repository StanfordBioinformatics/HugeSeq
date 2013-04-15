#!/bin/bash -eu

echo "*** Performing GFF Merging ***"

if [ $# -lt 2 ]
then
	echo "Usage: $0 <output> <GFF file to merge>..."
	exit 1
fi

output=`cd \`dirname $1\`; pwd`/`basename $1`
input=${output/.gff/}.raw.gff

shift

inputs=''
for i in $*
do
	inputs="$inputs `cd \`dirname $i\`; pwd`/`basename $i`"
done

cat $inputs > $input

features="`cut -f 3 $input | sort -u`"
sources="`cut -f 2 $input | sort -u`"

for feat in $features
do
	f=$input.${feat}
	grep "	$feat	" $input > $f
	for src in $sources
	do
		grep "	$src	" $f | mergeBed -i stdin > $f.$src
	done
	cat $f.* > $f
	> $f.dup
	> $f.uni
	intersectBed -a $f -b $f -f 0.5 -r -c | awk '{print $0 > ($NF>1? "'$f.dup'": "'$f.uni'")}'
	mergeBed -n -i $f.uni > $f.uni.merged
	mergeBed -n -i $f.dup > $f.dup.merged
	for i in $f.*.merged
	do
		dup="0"
		if [ -n "${i/*.uni.merged/}" ]; then dup="1"; fi
		awk -F '\t' '{qual="LowQual"; if ('$dup' && $4>=2) qual="PASS"; print $1"\tHugeSeq\t'$feat'\t"$2"\t"$3"\t"qual"\t.\t.\tEVENTS "$4}' $i
	done
	rm $f $f.*
done > $output
