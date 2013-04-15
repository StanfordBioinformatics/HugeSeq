#!/bin/bash -eu

echo "*** Recalibrating base quality ***"

if [ $# -lt 2 ]
then 
	echo "Usage: $0 <bam> <out>"
	exit 1
fi

f=`cd \`dirname $1\`; pwd`/`basename $1`
o=`cd \`dirname $2\`; pwd`/`basename $2`

echo ">>> Counting covariates"
java -Xms5g -Xmx5g -jar $GATK/GenomeAnalysisTK.jar \
	-T CountCovariates \
	-cov ReadGroupCovariate \
	-cov QualityScoreCovariate \
	-cov CycleCovariate \
	-cov DinucCovariate \
	-R $REF \
	--knownSites $SNP \
	-I $f \
	-recalFile ${o/.bam/.csv} \
	-et NO_ET \
        -K /srv/gs1/projects/snyder/cuiping/data/referencefiles/GATKkey/cuiping_stanford.edu.key

echo ">> Table recalibration"
if [ "`grep -v '#' ${o/.bam/.csv} | grep -v "EOF" | wc -l`" = "1" ]
then
	cp $f $o
else
	java -Xms5g -Xmx5g -jar $GATK/GenomeAnalysisTK.jar \
	-R $REF \
	-I $f \
	-o $o \
	-T TableRecalibration \
	-baq RECALCULATE \
	--doNotWriteOriginalQuals \
	-recalFile ${o/.bam/.csv} \
	-et NO_ET \
	-K /srv/gs1/projects/snyder/cuiping/data/referencefiles/GATKkey/cuiping_stanford.edu.key
fi

echo "*** Finished recalibrating base quality ***"
