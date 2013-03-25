#!/bin/env python

import sys, gzip
from util import *

if len(sys.argv) <= 2:
	print "Usage: %s <fasta/q> <reads per file>" % sys.argv[0]
	exit(1)


ifile = File(sys.argv[1])
reads = int(sys.argv[2])
lines = reads * 4

gz = True if ifile.ext == "gz" else False

input = gzip.open(ifile.path, "rb") if gz else open(ifile.path, "r")

c = 0
i = 0
output = None
for l in input:
	i+=1
	if i%lines==1:
		c += 1
		if l.startswith(">"):
			lines = reads * 2
		ofilename = ifile.absprefix if gz else ifile.path
		ofile = File((ofilename+".S%06d"%c)+(".gz" if gz else ""))
		if output is not None:
			output.flush()
			output.close()
                output = gzip.open(ofile.path, "wb") if gz else open(ofile.path, "w")
		
	output.write(l)

if output is not None:
	output.flush()
	output.close()
input.close()
