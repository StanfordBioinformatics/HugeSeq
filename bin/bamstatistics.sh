#!/bin/bash -e
#Usage: $0 <bam1> [bam2...]

awkopt='BEGIN {m=0;u=0;} {m+=$3; u+=$4} END {t=m+u; mp=0; if (m>0||u>0) mp=(m+0.0)/t; up=0; if (u>0) up=1-mp; x=t*101/3000000000; print m"\t"mp"\t"u"\t"up"\t"t"\t"x}'

s=0

for i in $*; 
do 
        let "s += 1"
	echo -e "header: chr-name\tchr-length\t#_mapped-reads\t#_unmapped-reads\n\n" >> ${s}.bam.stat;
	echo -e "${i/_?_pf.*/}\n`samtools idxstats $i`" >> ${s}.bam.stat;
	echo -e "\n" >> ${s}.bam.stat;
	echo -e "#_mapped-reads\t%_mapped-reads\t#_unmapped-reads\t%_unmapped-reads\t#_total-reads\ttotal-coverage" >> ${s}.bam.stat;
        echo -e "`awk \"$awkopt\" ${s}.bam.stat` \n" >> ${s}.bam.stat; 
done

echo -e "#_lanes\t#_mapped-reads\t%_mapped-reads\t#_unmapped-reads\t%_unmapped-reads\t#_total reads\ttotal-coverage" >> allbam.statistics
echo -e "Total $# lanes\t`(for i in $*; do samtools idxstats $i; done) | awk \"$awkopt\"`" >> allbam.statistics
