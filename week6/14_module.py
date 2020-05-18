#!/usr/bin/env python3

import gzip
import sys

# We have imported modules like math, sys, and gzip a few times
# You can also write and import your own modules

import biotools as bt

# The read_fasta() and gc() functions are in biotools.py (take a look)

for name, seq in bt.read_fasta('genome.fa.gz'):
	print(name, len(seq), bt.gc(seq))


