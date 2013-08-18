#!/bin/bash -eu

echo "*** Recalibrating base quality ***"

if [ $# -lt 2 ]
then 
	echo "Usage: $0 <bam> <out>"
	exit 1
fi

f=`cd \`dirname $1\`; pwd`/`basename $1`
o=`cd \`dirname $2\`; pwd`/`basename $2`

echo ">>> Counting covariates2"
java -Xms5g -Xmx5g -jar $GATK/GenomeAnalysisTK.jar \
	-T BaseRecalibrator \
	-I $f \
   	-R $REF \
	-o ${o/.bam/.grp} \
   	-knownSites $SNP \
	-cov ReadGroupCovariate \
	-cov QualityScoreCovariate \
	-cov CycleCovariate \
        -K /srv/gs1/software/gatk/GATKkey/stanford.edu.key
	#-cov DinucCovariate \
	#-recalFile ${o/.bam/.csv} \
	#-o ${o/.bam/.csv} \
   	#-o recal_data.table
	# output:? recal_data.grp
	#-et NO_ET \
#exit
#echo ">> Table recalibration"
#if [ "`grep -v '#' ${o/.bam/.csv} | grep -v "EOF" | wc -l`" = "1" ]
#then
#	cp $f $o
#else
	java -Xms5g -Xmx5g -jar $GATK/GenomeAnalysisTK.jar \
	-R $REF \
	-I $f \
	-o $o \
	-T PrintReads \
	-BQSR ${o/.bam/.grp} \
	-K /srv/gs1/software/gatk/GATKkey/stanford.edu.key
	#-T TableRecalibration \
	#-baq RECALCULATE \
	#--doNotWriteOriginalQuals \
	#-recalFile ${o/.bam/.csv} \
	#-et NO_ET \
#fi

echo "*** Finished recalibrating base quality ***"
