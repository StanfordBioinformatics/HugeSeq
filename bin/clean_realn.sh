#!/bin/bash -eu

echo "*** Realigning targeted regions ***"

if [ $# -lt 2 ]
then 
	echo "Usage: $0 <bam> <out>"
	exit 1
fi


f=`cd \`dirname $1\`; pwd`/`basename $1`
o=`cd \`dirname $2\`; pwd`/`basename $2`

optL=''
if [[ "$1" =~ ".*chr[^.]*\..*" ]]
then
	chr=`echo "$1" | sed 's/.*\(chr[^\.]*\)\..*/\1/'`
	optL="-L $chr"
fi

echo ">>> Determining (small) suspicious intervals which are likely in need of realignment"
java -Xms8g -Xmx8g -jar $GATK/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-I $f \
	-R $REF \
	-o ${o/.bam/.intervals} \
	-D $SNP $optL \
	-et NO_ET

echo ">>> Running the realigner over the targeted intervals"
java -Xms8g -Xmx8g -Djava.io.tmpdir=$TMP -jar $GATK/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-I $f \
	-R $REF \
	-o $o \
	-targetIntervals ${o/.bam/.intervals} \
	-D $SNP \
	-LOD 5 $optL \
	-et NO_ET

echo ">>> Fixing the mate pairs and order of the realigned reads"
java -Xms8g -Xmx8g -jar $PICARD/FixMateInformation.jar \
	TMP_DIR=$TMP \
	INPUT=$o \
	VALIDATION_STRINGENCY=SILENT \
	SORT_ORDER=coordinate

echo "*** Finished realigning targeted regions ***"
