#!/bin/bash

echo ">>> Writing reference calls"
echo "   " >> $LOGFILE
echo ">>> Writing reference calls" >> $LOGFILE

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

REFCALL_VCF=$PREFIX/$SAMPLE.refcalls.vcf

command="java -Xmx6g -Xms6g -jar $GATK/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R $REF \
   -V $START_VCF \
   -o $REFCALL_VCF \
   -selectType NO_VARIATION"

echo ">>> Select reference calls"
echo ">>> Select reference calls" >> $LOGFILE
echo ">>> $command &> $PREFIX/$SAMPLE.select.refcalls.log"
echo ">>> $command &> $PREFIX/$SAMPLE.select.refcalls.log" >> $LOGFILE
$command &> $PREFIX/$SAMPLE.select.refcalls.log 

echo ">>> Finished writine reference calls"
