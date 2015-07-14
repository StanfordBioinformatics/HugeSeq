#!/bin/bash -eu

echo "*** Sorting BAM by position ***"

if [ $# -lt 2 ]
then 
	echo "Usage: $0 <bam> <sorted> [memory in GB]"
	exit 1
fi

f=`cd \`dirname $1\`; pwd`/`basename $1`
o=`cd \`dirname $2\`; pwd`/`basename $2`

if [ "$2" = '-' ]
then
	o=$f.sorted
fi

gmem=5
if [ $# -gt 2 ]
then
	gmem=$3
fi

command="java -Xms${gmem}g -Xmx${gmem}g -jar $PICARD/SortSam.jar \
	TMP_DIR=$TMP \
	INPUT=$f \
	OUTPUT=$o \
	MAX_RECORDS_IN_RAM=$(($gmem*250000)) \
	VALIDATION_STRINGENCY=SILENT \
	SORT_ORDER=coordinate"

echo ">>> Sorting on BAM $f"
echo ">>> Sorting on BAM $f" >> $LOGFILE
echo ">>> $command"
echo ">>> $command" >> $LOGFILE
$command

if [ "$2" = '-' ]
then
	mv $o $f
	o=$f
fi

echo "*** Finished sorting BAM by position ***"
