#!/bin/bash

echo "Performing VQSR on indels" 
echo "   " >> $LOGFILE
echo "Performing VQSR on indels" >> $LOGFILE

set -e

if [ $# -lt 1 ]
then
        echo "Usage: $0 <vcf> vqsr targeted"
        exit 1
fi

START_VCF=`cd \`dirname $1\`; pwd`/`basename $1`
PREFIX=`dirname $START_VCF`
SUFFIX=`basename $START_VCF`
SAMPLE=${SUFFIX/.gatk.vcf/}

INDEL_VCF=$PREFIX/$SAMPLE.indel.vcf
INDEL_RECAL_VCF=$PREFIX/$SAMPLE.vqsr.indel.vcf
INDEL_RECAL=$PREFIX/$SAMPLE.tmp.indel.vcf
INDEL_TRANCHES=$PREFIX/$SAMPLE.tranches.gatk.indel.recal.csv
INDEL_RSCRIPT=$PREFIX/$SAMPLE.gatk.recal.indel.R

touch $INDEL_RECAL_VCF
doVQSR=$2
targeted=$3
    
if [ "$targeted" != "True" ]
then
	DP="-an DP"
fi

command="java -Xmx6g -Xms6g -jar $GATK/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R $REF \
   -V $START_VCF \
   -o $INDEL_VCF \
   -selectType INDEL"

echo ">>> Select variants for indel recalibration"
echo ">>> Select variants for indel recalibration" >> $LOGFILE
echo ">>> $command &> $PREFIX/$SAMPLE.select.indel.log"
echo ">>> $command > $PREFIX/$SAMPLE.select.indel.log" >> $LOGFILE
$command &> $PREFIX/$SAMPLE.select.indel.log

if [ "$doVQSR" != "True" ]
then
	command="cp $INDEL_VCF $INDEL_RECAL_VCF"
	echo ">>> Do not perform indel VQSR"
	echo ">>> Do not perform indel VQSR" >> $LOGFILE
	echo "$command > $PREFIX/$SAMPLE.vqsr.indel.log" >> $LOGFILE
	$command &> $PREFIX/$SAMPLE.vqsr.indel.log
	exit $?
fi

# http://gatkforums.broadinstitute.org/discussion/1259/what-vqsr-training-sets-arguments-should-i-use-for-my-specific-project
# -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.b37.vcf \
# -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \

command="java -Xmx6g -Xms6g -jar $GATK/GenomeAnalysisTK.jar \
    -T VariantRecalibrator \
    -R $REF \
    -input $INDEL_VCF \
    -resource:mills,known=true,training=true,truth=true,prior=12.0 $MILLS_1K_GOLD_INDELS \
    -an FS $DP \
    -an MQRankSum \
    -an ReadPosRankSum \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    -mode INDEL \
    -recalFile $INDEL_RECAL \
    -tranchesFile $INDEL_TRANCHES \
    --maxGaussians 4 \
    -rscriptFile $INDEL_RSCRIPT"

echo ">>> Train recalibration for indels"
echo ">>> Train recalibration for indels" >> $LOGFILE
echo ">>> $command &> $PREFIX/$SAMPLE.recalibrate.indel.log"
echo ">>> $command > $PREFIX/$SAMPLE.recalibrate.indel.log" >> $LOGFILE
$command &> $PREFIX/$SAMPLE.recalibrate.indel.log

# --ts_filter_level 99.0 \
# http://gatkforums.broadinstitute.org/discussion/1259/what-vqsr-training-sets-arguments-should-i-use-for-my-specific-project
# 99.9

command="java -Xmx6g -Xms6g -jar $GATK/GenomeAnalysisTK.jar \
   -T ApplyRecalibration \
   -R $REF \
   -input $INDEL_VCF \
   --ts_filter_level 99.0 \
   -tranchesFile $INDEL_TRANCHES \
   -recalFile $INDEL_RECAL \
   -o $INDEL_RECAL_VCF \
   --mode INDEL"

echo ">>> Apply recalibration for indels"
echo ">>> Apply recalibration for indels" >> $LOGFILE
echo ">>> $command &> $PREFIX/$SAMPLE.apply.indel.log"
echo ">>> $command > $PREFIX/$SAMPLE.apply.indel.log" >> $LOGFILE
$command &> $PREFIX/$SAMPLE.apply.indel.log

echo "Finished performing VQSR on indels" 
#rm $SNP_VCF $SNP_RECAL $SNP_TRANCHES $PREFIX/*.log $INDEL_VCF $INDEL_RECAL $INDEL_TRANCHES
