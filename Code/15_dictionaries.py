#!/usr/bin/env python3

import math
import random
import time
import gzip
import string

# Move the triple quotes downward to uncover each segment of code




# A list is constructed with square brackets
# and indexed by integers in square brackets

alist = ['cat', 'dog', 'fox']
print(alist[1])

# A dictionary is constructed with curly brackets
# Individual items have key:value pairs where the key is usually a string
# Dictionaries are usually indexed by text in square brackets

adict = {'cat':'meow', 'dog':'woof', 'fox':'????'}
print(adict['dog'])

# It is much faster to look up values in a dictionary than list

# List of 1 million elements
alist = []
for i in range(1000000): alist.append(0)
alist[999999] = 1
t0 = time.perf_counter()
c = 0
for i in range(100):
	if 1 in alist: c += 1
t1 = time.perf_counter()
et_list = t1 - t0
print('list', et_list)

# Dictionary of 1 million elements
adict = {}
for i in range(1000000): adict[i] = 0
adict[999999] = 1
t0 = time.perf_counter()
c = 0
for i in range(100):
	if 1 in adict: c += 1
t1 = time.perf_counter()
et_dict = t1 - t0
print('dict', et_dict)
print('speed up', et_list/et_dict)

# Dictionaries can be used for efficent look-up
# This is much faster than if-elif-else
gcode = {
	'AAA' : 'K',	'AAC' : 'N',	'AAG' : 'K',	'AAT' : 'N',
	'ACA' : 'T',	'ACC' : 'T',	'ACG' : 'T',	'ACT' : 'T',
	'AGA' : 'R',	'AGC' : 'S',	'AGG' : 'R',	'AGT' : 'S',
	'ATA' : 'I',	'ATC' : 'I',	'ATG' : 'M',	'ATT' : 'I',
	'CAA' : 'Q',	'CAC' : 'H',	'CAG' : 'Q',	'CAT' : 'H',
	'CCA' : 'P',	'CCC' : 'P',	'CCG' : 'P',	'CCT' : 'P',
	'CGA' : 'R',	'CGC' : 'R',	'CGG' : 'R',	'CGT' : 'R',
	'CTA' : 'L',	'CTC' : 'L',	'CTG' : 'L',	'CTT' : 'L',
	'GAA' : 'E',	'GAC' : 'D',	'GAG' : 'E',	'GAT' : 'D',
	'GCA' : 'A',	'GCC' : 'A',	'GCG' : 'A',	'GCT' : 'A',
	'GGA' : 'G',	'GGC' : 'G',	'GGG' : 'G',	'GGT' : 'G',
	'GTA' : 'V',	'GTC' : 'V',	'GTG' : 'V',	'GTT' : 'V',
	'TAA' : '*',	'TAC' : 'Y',	'TAG' : '*',	'TAT' : 'Y',
	'TCA' : 'S',	'TCC' : 'S',	'TCG' : 'S',	'TCT' : 'S',
	'TGA' : '*',	'TGC' : 'C',	'TGG' : 'W',	'TGT' : 'C',
	'TTA' : 'L',	'TTC' : 'F',	'TTG' : 'L',	'TTT' : 'F',
}
dna = 'ATGGTTCGCCAAGGAAACCTCCATCCCATGTNTCGGAGAATCTCTCCCGCTCCGACCGCT'
# note the N right here               ^
pro = []
for i in range(0, len(dna), 3):
	codon = dna[i:i+3]
	if codon in gcode: pro.append(gcode[codon])
	else:              pro.append('X') # how to deal with Ns and such
print(''.join(pro))

# Dictionaries can be used to count things, like k-mers
k = 2
count = {}
for i in range(len(dna) -k +1):
	kmer = dna[i:i+k]
	if kmer in count: count[kmer] += 1
	else:             count[kmer] = 1
for kmer in count:
	print(kmer, count[kmer])
