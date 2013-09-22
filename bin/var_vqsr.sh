#!/bin/bash

set -e

if [ $# -lt 1 ]
then
        echo "Usage: $0 <vcf>"
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

INDEL_VCF=$PREFIX/$SAMPLE.indel.vcf
INDEL_RECAL_VCF=$PREFIX/$SAMPLE.vqsr.indel.vcf
INDEL_RECAL=$PREFIX/$SAMPLE.tmp.indel.vcf
INDEL_TRANCHES=$PREFIX/$SAMPLE.tranches.gatk.indel.recal.csv
INDEL_RSCRIPT=$PREFIX/$SAMPLE.gatk.recal.indel.R

touch $SNP_RECAL_VCF
touch $INDEL_RECAL_VCF

java -Xmx6g -Xms6g -jar $GATK/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R $REF \
   -V $START_VCF \
   -o $SNP_VCF \
   -selectType SNP &> $PREFIX/$SAMPLE.select.snv.log

java -Xmx6g -Xms6g -jar $GATK/GenomeAnalysisTK.jar \
   -T VariantRecalibrator \
   -R $REF \
   -input $SNP_VCF \
   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $VQSR_HAPMAP \
   -resource:omni,known=false,training=true,truth=true,prior=12.0 $VQSR_OMNI \
   -resource:1000G,known=false,training=true,truth=false,prior=10.0 $VQSR_SNPS \
   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $SNP \
   -an DP \
   -an QD \
   -an FS \
   -an MQRankSum \
   -an ReadPosRankSum \
   -mode SNP \
   --percentBadVariants 0.05 \
   --maxGaussians 4 \
   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
   -recalFile $SNP_RECAL \
   -tranchesFile $SNP_TRANCHES \
   -rscriptFile $SNP_RSCRIPT &> $PREFIX/$SAMPLE.recalibrate.snv.log
   #--maxGaussians 4 \
   #--minNumBadVariants 1000
   #-mode SNP $VqsrMinNumBadVariants \

java -Xmx3g -Xms3g -jar $GATK/GenomeAnalysisTK.jar \
   -T ApplyRecalibration \
   -R $REF \
   -input $SNP_VCF \
   --ts_filter_level 99.0 \
   -tranchesFile $SNP_TRANCHES \
   -recalFile $SNP_RECAL \
   -o $SNP_RECAL_VCF \
   --mode SNP &> $PREFIX/$SAMPLE.apply.snv.log

java -Xmx6g -Xms6g -jar $GATK/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R $REF \
   -V $START_VCF \
   -o $INDEL_VCF \
   -selectType INDEL &> $PREFIX/$SAMPLE.select.indel.log

java -Xmx6g -Xms6g -jar $GATK/GenomeAnalysisTK.jar \
    -T VariantRecalibrator \
    -R $REF \
    -input $INDEL_VCF \
    -resource:mills,known=true,training=true,truth=true,prior=12.0 $VQSR_INDELS \
    -an DP \
    -an FS \
    -mode INDEL $VqsrMinNumBadVariants \
    -percentBad 0.01 \
    -an ReadPosRankSum \
    -an MQRankSum \
    --percentBadVariants 0.05 \
    --maxGaussians 4 \
    -minNumBad 1000 \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    -recalFile $INDEL_RECAL \
    -tranchesFile $INDEL_TRANCHES \
    -rscriptFile $INDEL_RSCRIPT &> $PREFIX/$SAMPLE.recalibrate.indel.log 

java -Xmx6g -Xms6g -jar $GATK/GenomeAnalysisTK.jar \
   -T ApplyRecalibration \
   -R $REF \
   -input $INDEL_VCF \
   --ts_filter_level 99.0 \
   -tranchesFile $INDEL_TRANCHES \
   -recalFile $INDEL_RECAL \
   -o $INDEL_RECAL_VCF \
   --mode INDEL &> $PREFIX/$SAMPLE.apply.indel.log 

#rm $SNP_VCF $SNP_RECAL $SNP_TRANCHES $PREFIX/*.log $INDEL_VCF $INDEL_RECAL $INDEL_TRANCHES
