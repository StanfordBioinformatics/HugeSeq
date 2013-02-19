#!/bin/bash -eu

echo "*** Removing duplicates ***"

if [ $# -lt 2 ]
then 
	echo "Usage: $0 <bam> <out> [remove, default: true]"
	exit 1
fi

rmdup="true"
if [ $# -gt 2 ]
then
	rmdup=$3
fi

f=`cd \`dirname $1\`; pwd`/`basename $1`
o=`cd \`dirname $2\`; pwd`/`basename $2`

echo ">>> Marking duplicates"
java -Xms5g -Xmx5g -jar $PICARD/MarkDuplicates.jar \
	TMP_DIR=$TMP \
	I=${f} \
	O=${o} \
	M=${o/.bam/.metrics} \
	VALIDATION_STRINGENCY=SILENT \
	ASSUME_SORTED=true \
	REMOVE_DUPLICATES=$rmdup

echo "*** Finished removing duplicates ***"
