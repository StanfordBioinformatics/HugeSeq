#!/bin/env python

## Fixes BAM index file name to match what GATK expects.

import os
import sys


orig_filename = sys.argv[1]

if not os.path.exists(orig_filename):
    print "Cannot find file", sys.argv[1]
    raise SystemExit(1)
if not orig_filename.endswith('.bam.bai'):
    print "Filename doesn't end with .bam.bai"
    raise SystemExit(0)

new_filename = orig_filename[:-8] + '.bai'

if os.path.exists(new_filename):
    print "Existing file %s renamed to %s.old" % (new_filename, new_filename)
    os.rename(new_filename, new_filename + '.old')

#os.rename(orig_filename, new_filename)
os.link(orig_filename, new_filename)
