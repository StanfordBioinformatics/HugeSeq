#!/bin/bash -eu

echo "*** Marking duplicates ***"

if [ $# -lt 2 ]
then 
	echo "Usage: $0 <bam> <out> [remove, default: false]"
	exit 1
fi

input=`cd \`dirname $1\`; pwd`/`basename $1`
output=`cd \`dirname $2\`; pwd`/`basename $2`

cp $input $output

command="java -Xms5g -Xmx5g -jar $PICARD/MarkDuplicates.jar \
	TMP_DIR=${TMP} \
	I=${input} \
	O=${output} \
	M=${output/.bam/.metrics} \
	VALIDATION_STRINGENCY=SILENT \
	ASSUME_SORTED=true \
	REMOVE_DUPLICATES=false"

echo ">>> Marking duplicates"
echo ">>> Marking duplicates" >> $LOGFILE
echo ">>> $command"
echo ">>> $command" >> $LOGFILE
$command

echo "*** Finished marking duplicates ***"
