#!/bin/sh

if [ $# -lt 2 ]
then
	echo "Usage: $0 <input fasta/q> <link to input>"
	exit 1
fi

i=`cd \`dirname $1\`; pwd`/`basename $1`
l=`cd \`dirname $2\`; pwd`/`basename $2`

if [ ! -e $l -a $i != $l ]
then
	echo ">> Creating link to input sequence file"
	echo "-- Input: $i"
	echo "-- Link : $l"
	ln -s $i $l

	if [[ $i =~ ".*\.bam" ]]
	then
		echo "-- Link : $l.bai"
		ln -s $i.bai $l.bai
	fi
fi
