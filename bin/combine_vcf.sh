#!/bin/bash -eu

if [ $# -lt 2 ]
then 
	echo "Usage: $0 <bam> <out>"
	exit 1
fi

out=`cd \`dirname $1\`; pwd`/`basename $1`
shift

noop=$1
shift

zip=$1
shift

inputs=`cd \`dirname $1\`; pwd`/`basename $1`
shift
for i in $*
do
       	inputs="$inputs --variant `cd \`dirname $i\`; pwd`/`basename $i`"
done

if [[ "$noop" == "False" ]]
then
	command="java -Xmx8g -Xms8g -jar $GATK/GenomeAnalysisTK.jar \
	   -R $REF \
	   -T CombineVariants \
	   --variant $inputs \
	   -o $out"

	echo ">>> Combining VCFs"
	echo ">>> Combining VCFs" >> $LOGFILE
	echo ">>> $command"
	echo ">>> $command" >> $LOGFILE
	$command

	if [[ "$zip" == "True" ]]
	then
		command="bgzip -fc $out"
	        echo ">>> Zipping VCF"
        	echo ">>> $command"
	        echo ">>> $command" >> $LOGFILE
        	$command &> $out.gz

	        command="tabix -p vcf $out.gz"
        	echo ">>> Indexing VCF using tabix"
	        #echo ">>> Indexing VCF using tabix" >>> $LOGFILE
        	echo ">>> $command"
	        echo ">>> $command" >> $LOGFILE
        	$command
	fi
fi
echo "*** Finished combining VCFs ***"
