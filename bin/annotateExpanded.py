#!/bin/env python

import sys
import os
import re

if len(sys.argv) <= 2:
	print 'usage: <output file> <vcf or gff file>'
	exit(1)
	
path = os.environ['ANNOVAR']
#path = /srv/gs1/projects/snyder/cuiping/app/annovar/
output = sys.argv[1]
input = sys.argv[2]

isVCF=False
if input.endswith('.vcf') or input.endswith('.vcf.gz'):
	isVCF=True
	vcf = input
	avinput = vcf + '.avinput'
	avoutput = vcf + '.avoutput'
	os.system('less %s | %s/convert2annovar.pl -format vcf4 - > %s' % (vcf, path, avinput))
elif input.endswith('.gff'):
	gff = input
	avinput = gff + '.avinput'
	avoutput = gff + '.avoutput'
	gff_file = open(gff, 'read')
	out_file = open(avinput, 'w')
	for line in gff_file:
		gffCols = line.split('\t')
		out_file.write(gffCols[0]+'\t'+gffCols[3]+'\t'+gffCols[4]+'\t0\t0\n')
	out_file.flush()
	out_file.close()
else:
	print >> sys.stderr, "Unknown input format: "+input
	exit(1)

# use user-define output name for now
avoutput=output

print 'Annotating variants with hg19 UCSC knownGene...\n'
#os.system('%s/annotate_variation.pl --geneanno --buildver hg19 -dbtype knownGene --separate %s %s/humandb/' %(path, avinput, path))
/usr/bin/perl $path/annotate_variation.pl --geneanno --buildver hg19 -dbtype knownGene --separate $avinput $path/humandb/

print 'Annotating variants with hg19 RMSK...\n'
#os.system('%s/annotate_variation.pl -regionanno --buildver hg19 -dbtype gff3 -gff3dbfile hg19_rmsk.gff %s %s/humandb/' %(path, avinput, path))
/usr/bin/perl $path/annotate_variation.pl -regionanno --buildver hg19 -dbtype gff3 -gff3dbfile hg19_rmsk.gff $avinput $path/humandb/

if isVCF:
	print 'Annotating variants with sift scores using hg19...\n'
#	os.system('%s/annotate_variation.pl --filter --sift_threshold 0 --buildver hg19 --separate -dbtype avsift %s %s/humandb/' %(path, avinput, path))
        /usr/bin/perl $path/annotate_variation.pl --filter --sift_threshold 0 --buildver hg19 --separate -dbtype avsift $avinput $path/humandb/

	print 'Annotating variants with dbSNP135...\n'
#	os.system('%s/annotate_variation.pl --filter --buildver hg19 --dbtype snp135 %s %s/humandb/' %(path, avinput, path))
        /usr/bin/perl $path/annotate_variation.p --filter --buildver hg19 --dbtype snp135 $avinput $path/humandb/

	#print 'Annotating variants with wgEncodeBroadHistoneGm12878H3k4me1StdPk...\n'
	#os.system('%s/annotate_variation.pl -regionanno --buildver hg19 -dbtype wgEncodeBroadHistoneGm12878H3k4me1StdPk -scorecolumn 5 %s %s/humandb/' %(path, avinput, path))

        print 'Annotating variants with 1000g...\n'
	/usr/bin/perl $path/annotate_variation.pl --filter --buildver hg19 --dbtype 1000g2012feb_all $avinput $path/humandb/


exonic_file = avinput + '.exonic_variant_function'
function_file = avinput + '.variant_function'
sift_file = avinput + '.hg19_avsift_dropped'
dbsnp_file = avinput + '.hg19_snp132_dropped'
#meth_file = avinput + '.hg19_wgEncodeBroadHistoneGm12878H3k4me1StdPk'
rmsk_file = avinput + ".hg19_gff3"
1000g_file = avinput + ".hg19_ALL.sites.2012_02_dropped"

def makeDict(filename, chrCol, startCol, endCol, valueCols, isGFF=False):
	file = open(filename, 'r')
	dic = {}
	for line in file:
		cols = line.split('\t')
		key=(cols[chrCol], cols[startCol], cols[endCol])
		if key not in dic:
			dic[key] = []
			for valueCol in valueCols:
				dic[key].append([])
		for i in range(len(valueCols)):
			value=cols[valueCols[i]]
			if isGFF:
				value=value.split(";")[-1].split("=")[-1]
			dic[key][i].append(value)
	return dic

function = makeDict(function_file, 2, 3, 4, [1, 0])
rmsk = makeDict(rmsk_file, 2, 3, 4, [1], True)
if isVCF:
	exonic = makeDict(exonic_file, 3, 4, 5, [1,2])
	dbsnp = makeDict(dbsnp_file, 2, 3, 4, [1])
	#meth = makeDict(meth_file, 2, 3, 4, [1])
	sift = makeDict(sift_file, 2, 3, 4, [1])

AVINPUT = open(avinput, 'r')
AVOUTPUT = open(avoutput, 'w')
AVOUTPUT.write("#chr\tstart\tend\tgene_name\ttype\trmsk");
if isVCF:
	AVOUTPUT.write("\tSIFT\tconsequence\tmutation_info\tdbsnp132")

AVOUTPUT.write("\n")

def write(dic, key, nvalueCols=1):
	if key in dic: 
		for value in dic[key]:
			AVOUTPUT.write('\t'+";".join(value))
	else:		
		AVOUTPUT.write('\t.'*nvalueCols)

for line in AVINPUT:
	if re.match('([0-9A-Za-z]+)\s+(\d+)', line):
		splitline = line.split('\t')
		key=(splitline[0], splitline[1],splitline[2])
		AVOUTPUT.write(splitline[0]+'\t'+splitline[1]+'\t'+splitline[2])
		write(function, key, 2)
		write(rmsk, key)
		if isVCF:
			write(sift, key)
			write(exonic, key, 2)
			#write(meth, key)
			write(dbsnp, key)
		AVOUTPUT.write('\n')

AVINPUT.flush();
AVOUTPUT.flush();

AVINPUT.close();
AVOUTPUT.close();

#`rm $avinput.*`
