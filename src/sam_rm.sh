#!/bin/sh

echo "*** SAM/BAM Removel ***"

for i in $*
do
	if [ -z "${i/*.bam/}" -o -z "${i/*.sam/}" ]
	then
		echo ">> Removing SAM/BAM file: $i"
		rm -f $i $i.bai
	fi
done

echo "*** Finished SAM/BAM Removal ***"
