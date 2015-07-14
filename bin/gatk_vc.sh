#!/bin/bash -eu

echo "*** Performing SNV discovery and genotyipng using GATK ***"
echo "   " >> $LOGFILE
echo "*** Performing SNV discovery and genotyipng using GATK ***" >> $LOGFILE

if [ $# -lt 4 ]
then 
	echo "Usage: $0 output reference_calls snp_hapcaller indel_hapcaller bam ..."
	exit 1
fi

o=`cd \`dirname $1\`; pwd`/`basename $1`
shift
output_hc=$o.hc.vcf
output_gc=$o.gc.vcf

snp_vcf=$o.snp
indel_vcf=$o.indel

capture=$1
shift
echo $capture

reference_calls=$1
shift

snp_hap=$1
shift
indel_hap=$1
shift

optL=''
if [ $# -eq 1 ]
then
	echo $1
        if [[ "$1" =~ .*chr[^\.]*\..* ]]
        then
                chr=`echo "$1" | sed 's/.*\(chr[^\.]*\)\..*/\1/'`
                optL="-L $chr"
        fi
fi

f=''
for i in $*
do
        f="$f -I `cd \`dirname $i\`; pwd`/`basename $i`"
done

log=${o/gatk.vcf/}vc.log

run_hc="False"
run_gc="False"
if [[ "$snp_hap" == "True" ]] || [[ "$indel_hap" == "True" ]] 
then
	run_hc="True"
fi

if [[ "$snp_hap" != "True" ]] || [[ "$indel_hap" != "True" ]] 
then
	run_gc="True"
fi

NO_VARIATION=""
OUTPUT_MODE_UG=""
OUTPUT_MODE_HC=""
if [[ "$reference_calls" == "True" ]]
then
	OUTPUT_MODE_UG="--output_mode emit_all_sites"
	OUTPUT_MODE_HC="-ERC BP_RESOLUTION"
	NO_VARIATION="-selecttype no_variation"
fi

CAPTURE=""
if [[ "$capture" != "False" ]]
then
	#CAPTURE="-L $capture"
	captures=$(echo $capture | tr "[" "\n")
	captures=$(echo $captures | tr "]" "\n")
	captures=$(echo $captures | tr "," "\n")
	#echo "here" $captures
	
	for capt in $captures
	do
		#echo "HELLO" $capt
		if [[ "$CAPTURE" != "" ]]
		then 
			CAPTURE="$CAPTURE -L $capt"
		else
			CAPTURE="-L $capt"
		fi
	done
fi

gc_command="java -Xmx8g -Xms8g -Djava.io.tmpdir=$TMP -jar $GATK/GenomeAnalysisTK.jar \
   -T UnifiedGenotyper \
   -R $REF \
   $f \
   --dbsnp $DBSNP \
   -o $output_gc $optL $CAPTURE \
   -stand_call_conf 20.0 \
   -stand_emit_conf 10.0 \
   -gt_mode DISCOVERY $OUTPUT_MODE_UG \
   --genotype_likelihoods_model BOTH \
   -nct 4"

hc_command="java -Xmx12g -Xms12g -Djava.io.tmpdir=$TMP -jar $GATK/GenomeAnalysisTK.jar \
     -T HaplotypeCaller \
     -R $REF \
     $f \
     --dbsnp $DBSNP \
     -o $output_hc $optL $CAPTURE \
     -stand_call_conf 20.0 \
     -stand_emit_conf 10.0 \
     --genotyping_mode DISCOVERY $OUTPUT_MODE_HC \
     -nct 4"

# https://gatkforums.broadinstitute.org/discussion/3115/emit-all-sites-in-haplotypecaller
# http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html

if [[ "$run_gc" == "True" ]]
then
	echo ">>> Running the GATK's UnifiedGenotyper"
	echo ">>> Running the GATK's UnifiedGenotyper" >> $LOGFILE
	echo ">>> $gc_command"
	echo ">>> $gc_command" >> $LOGFILE
	$gc_command &>> $log 
fi

if [[ "$run_hc" == "True" ]]
then
	echo ">>> Running the GATK's Haplotyper"
	echo ">>> Running the GATK's Haplotyper" >> $LOGFILE
	echo ">>> $hc_command"
	echo ">>> $hc_command" >> $LOGFILE
	$hc_command &>> $log
fi

use_gc="False"
use_hc="False"

if [[ "$snp_hap" == "True" ]] && [[ "$indel_hap" == "True" ]]
then
	# select both snps and indels from HC
	select_indel_from=$output_hc
	select_snp_from=$output_hc
	use_hc="True"
fi

if [[ "$snp_hap" == "True" ]] && [[ "$indel_hap" != "True" ]]
then
	# select snps from HC and indels from GC
	select_indel_from=$output_gc
	select_snp_from=$output_hc
fi

if [[ "$snp_hap" != "True" ]] && [[ "$indel_hap" == "True" ]]
then
	# select snps from GC and indels from HC
	select_indel_from=$output_hc
	select_snp_from=$output_gc
fi

if [[ "$snp_hap" != "True" ]] && [[ "$indel_hap" != "True" ]]
then
	# select snps and indels from GC
	select_indel_from=$output_gc
	select_snp_from=$output_gc
	use_gc="True"
fi

if [[ "$use_gc" == "True" ]]
then
	command="mv $output_gc $o"
	echo ">>> Choosing both SNP and INDEL from UG output: $select_snp_from"
	echo ">>> Choosing both SNP and INDEL from UG output: $select_snp_from" >> $LOGFILE
	echo ">>> $command"
	echo ">>> $command" >> $LOGFILE
	$command &>> $log

elif [[ "$use_hc" == "True" ]]
then
	command="mv $output_hc $o"
	echo ">>> Choosing both SNP and INDEL from HC output: $select_snp_from"
	echo ">>> Choosing both SNP and INDEL from HC output: $select_snp_from" >> $LOGFILE
	echo ">>> $command"
	echo ">>> $command" >> $LOGFILE
	$command &>> $log
else

	command="java -Xmx8g -Xms8g -jar $GATK/GenomeAnalysisTK.jar \
		-T SelectVariants \
		-R $REF \
		-V $select_snp_from \
		-o $snp_vcf \
		-selectType MNP \
		-selectType SNP $NO_VARIANTION"
	
		echo ">>> Selecting SNPs from $select_snp_from"
		echo ">>> Selecting SNPs from $select_snp_from" >> $LOGFILE
		echo ">>> $command"
		echo ">>> $command" >> $LOGFILE
		$command &>> $log

	command="java -Xmx8g -Xms8g -jar $GATK/GenomeAnalysisTK.jar \
		-T SelectVariants \
		-R $REF \
		-V $select_indel_from \
		-o $indel_vcf \
		-selectType INDEL"

		echo ">>> Selecting Indels from $select_indel_from"
		echo ">>> Selecting Indels from $select_indel_from" >> $LOGFILE
		echo ">>> $command"
		echo ">>> $command" >> $LOGFILE
		$command &>> $log

	command="java -Xmx8g -Xms8g -jar $GATK/GenomeAnalysisTK.jar \
		-R $REF \
		-T CombineVariants \
		--variant $snp_vcf \
		--variant $indel_vcf \
		-o $o"

		echo ">>> Combining SNP and Indel VCFs"
		echo ">>> Combining SNP and Indel VCFs" >> $LOGFILE
		echo ">>> $command"
		echo ">>> $command" >> $LOGFILE
		$command &>> $log
fi
	
echo "*** Finished SNV discovery and genotyipng using GATK ***"

