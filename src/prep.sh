#!/bin/sh

#export >> /srv/gs1/projects/snyder/cuiping/data/test/xxx


if [ $# -lt 2 ]
then
	echo "Usage: $0 <input fasta/q> <link to input>"
	exit 1
fi

i=`cd \`dirname $1\`; pwd`/`basename $1`
l=`cd \`dirname $2\`; pwd`/`basename $2`

echo ">> Creating link to input sequence file"
echo "-- Input: $i"
echo "-- Link : $l"

if [ ! -e $l ]
then
	ln -s $i $l
fi
