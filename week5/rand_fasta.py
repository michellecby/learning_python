#!/usr/bin/env python3

import gzip
import sys
import math
import random

# Write a program that finds creates random fasta files
# Create a function that makes random DNA sequences
# Parameters include length and frequencies for A, C, G, T
# Command line: python3 rand_fasta.py <count> <min> <max> <a> <c> <g> <t>

assert(len(sys.argv) == 8)
count = int(sys.argv[1])
min = int(sys.argv[2])
max = int(sys.argv[3])
a_ratio = float(sys.argv[4])
c_ratio = float(sys.argv[5])
g_ratio = float(sys.argv[6])
t_ratio = float(sys.argv[7])

assert(count >= 1)
assert(min >= 1)
assert(max >= min)
assert(a_ratio >= 0)
assert(c_ratio >= 0)
assert(t_ratio >= 0)
assert(g_ratio >= 0)
assert(math.isclose((a_ratio + c_ratio + g_ratio + t_ratio), 1))



def random_seq(values):
	r = random.randint(min, max)
	sequence = []
	for nt in range(r):
		r_nt = random.random()
		if r_nt <= a_ratio:
			sequence.append('A')
		elif r_nt <= c_ratio + a_ratio:
			sequence.append('C')
		elif r_nt <= g_ratio + c_ratio + a_ratio:
			sequence.append('G')
		else:
			sequence.append('T')
	return ''.join(sequence)

requires = [min, max, a_ratio, c_ratio, g_ratio, t_ratio]

for seq in range(count):
	print(f'>seq{seq}')
	print(random_seq(requires))

#print(random_seq(requires))
"""
python3 rand_fasta.py 5 2 10 0.25 0.25 0.25 0.25
python3 rand_fasta.py 15 200 1000 0.35 0.15 0.15 0.35  > output.fasta
"""
