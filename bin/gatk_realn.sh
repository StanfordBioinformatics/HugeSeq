#!/bin/bash -eu

echo "*** Realigning targeted regions ***"
echo "   " >> $LOGFILE
echo "*** Realigning targeted regions ***" >> $LOGFILE

if [ $# -lt 2 ]
then 
	echo "Usage: $0 <bam> <out>"
	exit 1
fi

f=`cd \`dirname $1\`; pwd`/`basename $1`
o=`cd \`dirname $2\`; pwd`/`basename $2`
relax_realign=$3

optL=''
if [[ "$1" =~ .*chr[^\.]*\..* ]]
then
	chr=`echo "$1" | sed 's/.*\(chr[^\.]*\)\..*/\1/'`
	optL="-L $chr"
fi

RELAX=""
if [[ "$relax_realign" == "True" ]]
then
        RELAX="--defaultBaseQualities 0 --filter_bases_not_stored --filter_mismatching_base_and_quals --filter_reads_with_N_cigar"
fi

command="java -Xms8g -Xmx8g -jar $GATK/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-I $f \
	-R $REF \
	-o ${o/.bam/.intervals} $optL \
        -known $MILLS_1K_GOLD_INDELS \
	-known $GOLD_1K_INDELS \
	-et NO_ET \
        -K /srv/gs1/software/gatk/GATKkey/stanford.edu.key"
 
echo ">>> Determining (small) suspicious intervals which are likely in need of realignment"
echo ">>> Determining (small) suspicious intervals which are likely in need of realignment" >> $LOGFILE
echo ">>> $command"
echo ">>> $command" >> $LOGFILE
$command

command="java -Xms8g -Xmx8g -Djava.io.tmpdir=$TMP -jar $GATK/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-I $f \
	-R $REF \
	-o $o $optL \
	-targetIntervals ${o/.bam/.intervals} \
        -known $MILLS_1K_GOLD_INDELS \
        -known $GOLD_1K_INDELS \
	-et NO_ET $RELAX\
	-K /srv/gs1/software/gatk/GATKkey/stanford.edu.key"
 
echo ">>> Running the realigner over the targeted intervals"
echo ">>> Running the realigner over the targeted intervals" >> $LOGFILE
echo ">>> $command"
echo ">>> $command" >> $LOGFILE
$command

command="java -Xms8g -Xmx8g -jar $PICARD/FixMateInformation.jar \
	TMP_DIR=$TMP \
	INPUT=$o \
	VALIDATION_STRINGENCY=SILENT \
	SORT_ORDER=coordinate"

echo ">>> Fixing the mate pairs and order of the realigned reads"
echo ">>> Fixing the mate pairs and order of the realigned reads" >> $LOGFILE
echo ">>> $command"
echo ">>> $command" >> $LOGFILE
$command

echo "*** Finished realigning targeted regions ***"
