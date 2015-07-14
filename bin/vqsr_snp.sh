#!/bin/bash

echo ">>> Performing VQSR on SNPs"
echo "   " >> $LOGFILE
echo ">>> Performing VQSR on SNPs" >> $LOGFILE

set -e

if [ $# -lt 1 ]
then
        echo "Usage: $0 <vcf> vqsr"
        exit 1
fi

START_VCF=`cd \`dirname $1\`; pwd`/`basename $1`
PREFIX=`dirname $START_VCF`
SUFFIX=`basename $START_VCF`
SAMPLE=${SUFFIX/.gatk.vcf/}

SNP_VCF=$PREFIX/$SAMPLE.snp.vcf
SNP_RECAL_VCF=$PREFIX/$SAMPLE.vqsr.snp.vcf
SNP_RECAL=$PREFIX/$SAMPLE.tmp.snp.vcf
SNP_TRANCHES=$PREFIX/$SAMPLE.tranches.gatk.snp.recal.csv
SNP_RSCRIPT=$PREFIX/$SAMPLE.gatk.recal.snp.R

touch $SNP_RECAL_VCF
doVQSR=$2
targeted=$3
useHaplotypeScore=$4
noop=$5

if [ "$noop" == "True" ]
then
	echo "bye bye"
        exit $?
fi

HAPLOTYPSCORE=""
if [ "$useHaplotypeScore" == "True" ]
then
        HAPLOTYPSCORE="-an HaplotypeScore"
fi

DP="-an DP"
if [ "$targeted" == "True" ]
then
        DP="--maxGaussians 4"
fi

if [ "$doVQSR" != "True" ]
then
        command="cp $SNP_VCF $SNP_RECAL_VCF"
        echo ">>> Do not perform SNP VQSR"
        echo ">>> Do not perform SNP VQSR" >> $LOGFILE
        $command &> $PREFIX/$SAMPLE.vqsr.snp.log
        echo "$command &> $PREFIX/$SAMPLE.vqsr.snp.log" >> $LOGFILE
        exit $?
fi

command="java -Xmx6g -Xms6g -jar $GATK/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R $REF \
   -V $START_VCF \
   -o $SNP_VCF \
   -selectType SNP"

echo ">>> Select variants for SNP recalibration"
echo ">>> Select variants for SNP recalibration" >> $LOGFILE
echo ">>> $command &> $PREFIX/$SAMPLE.select.snp.log"
echo ">>> $command &> $PREFIX/$SAMPLE.select.snp.log" >> $LOGFILE
$command &> $PREFIX/$SAMPLE.select.snp.log 

command="java -Xmx6g -Xms6g -jar $GATK/GenomeAnalysisTK.jar \
   -T VariantRecalibrator \
   -R $REF \
   -input $SNP_VCF \
   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
   -resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI_1K \
   -resource:1000G,known=false,training=true,truth=false,prior=10.0 $GOLD_1K_SNPS \
   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
   -an QD \
   -an FS $DP \
   -an MQRankSum \
   -an ReadPosRankSum $HAPLOTYPESCORE \
   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
   -mode SNP \
   -recalFile $SNP_RECAL \
   -tranchesFile $SNP_TRANCHES \
   -nt 4 \
   --minNumBadVariants 5000 \
   -rscriptFile $SNP_RSCRIPT"

echo ">>> Train recalibrator for SNPs"
echo ">>> Train recalibrator for SNPs" >> $LOGFILE
echo ">>> $command &> $PREFIX/$SAMPLE.recalibrate.snp.log"
echo ">>> $command > $PREFIX/$SAMPLE.recalibrate.snp.log" >> $LOGFILE
$command &> $PREFIX/$SAMPLE.recalibrate.snp.log

command="java -Xmx3g -Xms3g -jar $GATK/GenomeAnalysisTK.jar \
   -T ApplyRecalibration \
   -R $REF \
   -input $SNP_VCF \
   --ts_filter_level 99.0 \
   -tranchesFile $SNP_TRANCHES \
   -recalFile $SNP_RECAL \
   -o $SNP_RECAL_VCF \
   --mode SNP"

echo ">>> Apply recalibration to SNPs"
echo ">>> Apply recalibration to SNPs" >> $LOGFILE
echo ">>> $command &> $PREFIX/$SAMPLE.apply.snp.log"
echo ">>> $command > $PREFIX/$SAMPLE.apply.snp.log" >> $LOGFILE
$command &> $PREFIX/$SAMPLE.apply.snp.log

echo ">>> Finished performing VQSR on SNPs"
