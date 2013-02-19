#!/bin/sh

if [ $# -lt 2 ]
then
	echo "Usage: $0 <output file> <BAM> ..."
	exit 1
fi

#split_script=`cd \`dirname $0\`; pwd`/vcf_split.py

o=`cd \`dirname $1\`; pwd`/`basename $1`
shift
f=''
for i in $*
do
        f="$f `cd \`dirname $i\`; pwd`/`basename $i`"
done

echo ">> Pile up the reads and call SNPs for $f"
samtools mpileup -uf $REF $f | bcftools view -vcg - > $o

echo ">> Apply filter to the VCF and output to $o"
t=${o/.vcf/}.filtered.vcf 
vcfutils.pl varFilter -D4000  $o > $t 
mv $t $o

