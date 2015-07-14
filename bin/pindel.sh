#!/bin/bash -e

echo "*** Calling SV using Split-Read Analysis: $PINDEL ***"

if [ $# -lt 2 ]
then
	echo "Usage: $0 <output> <bam>..."
	exit 1
fi


o=`cd \`dirname $1\`; pwd`/`basename $1`
odir=`dirname $1`
shift
p=${o/.gff/}

rpm_output=${p/.pindel/}.breakdancer.txt
rpm_cfg=${p/.pindel/}.breakdancer.cfg

chr="ALL"
if [ $# -eq 1 ]
then
	if [[ "$1" =~ ".*chr[^.]*\..*" ]]
	then
		chr=`echo "$1" | sed 's/.*\(chr[^\.]*\)\..*/\1/'`
	fi
fi

bOpt=''
if [ -e $rpm_output ]
then
	bOpt="-b $rpm_output"
fi

cd $odir
pin_cfg=$p.cfg
for bam in $*
do
	isize=''
	if [ -e $rpm_cfg ]
	then
		cfg=`grep $bam $rpm_cfg | head -n 1`
		isize=`echo "$cfg" | sed 's/.*mean:\([0-9]*\).*/\1/'`
	else
		isize=`samtools view -H $bam 2> /dev/null | grep "@RG" | grep "PI" | head -n 1 | sed 's/.*PI:\([0-9]*\).*/\1/'` 
	fi
	if [ -z "$isize" ]
	then
		isize=300
	fi

	sample=`samtools view -H $bam 2> /dev/null | grep "@RG" | grep "SM" | head -n 1| sed 's/.*SM:\([^\t]*\).*/\1/'`
	if [ -z "$sample" ]
	then
		sample="SAMPLE"
	fi

	echo -e "$bam\t$isize\t$sample"
done > $pin_cfg
pindel="$PINDEL/pindel -f $REF -i $pin_cfg -o $p -c $chr $bOpt"
echo ">> Running Pindel on $chr"
$pindel

echo ">> Converting output to GFF"

minsize=50

AWKOPT='{feat="Unknown"; if ($2=="D" || $2=="DI") feat="Deletion"; else if ($2=="I" || $2=="LI") feat="Insertion"; else if ($2=="INV") feat="Inversion"; else if ($2=="TD") feat="TandemDup"; start=$7; end=$8; if (feat!="Insertion") {start++; end--;} if ($3>='$minsize') print $5"\tPindel\t"feat"\t"start"\t"end"\t"$24"\t.\t.\tSize "$3"; nr.unique reads: (+"$17",-"$20"); ComScore: "int(sqrt(($17+1)*($20+1)*$24))}'

#echo -e "#Chr\tProgram\tSV-type\t\tstart\tend\tscore\tstrand\tframe\tattributes" > $o
if [ -e "$o" ]
then
        rm $o
fi

for po in ${p}_[^LBT]*
do
	grep "ChrID" $po | sed 's/D \([0-9]*\)	I \([0-9]*\)/DI \1:\2/' | sed 's/SUM_MS \([0-9]*\)	.*/SUM_MS \1/' | sed 's/NT[^	]*//' | awk "$AWKOPT"
done | sort -k1,1 -k4n -k5n -u >> $o

#echo ">> Archiving raw output"
#tar --remove-files -zcvf $p.tgz ${p}_*

echo "*** Finished calling SV using Split-Read Analysis: $PINDEL ***"
