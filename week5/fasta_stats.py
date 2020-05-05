#!/usr/bin/env python3

import gzip
import sys

# Write a program that computes typical stats for sequence files
# Command line:
#	python3 fasta_stats.py transcripts.fasta.gz
# Output:
#	See below

def read_fasta(filename):
	name = None
	seqs = []

	fp = None
	if filename == '-':
		fp = sys.stdin
	elif filename.endswith('.gz'):
		fp = gzip.open(filename, 'rt')
	else:
		fp = open(filename)

	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield(name, seq)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield(name, ''.join(seqs))
	fp.close()

protein_file = sys.argv[1]

sizes = []
a, c, g, t = 0, 0, 0, 0
for name, seq in read_fasta(protein_file):
	sizes.append(len(seq))
	for nt in seq:
		if nt == 'A':   a += 1
		elif nt == 'C': c += 1
		elif nt == 'G': g += 1
		elif nt == 'T': t += 1



sizes.sort()
min = sizes[0]
max = sizes[-1]
total_len = a + c + g + t
a_ratio = a/total_len
c_ratio = c/total_len
g_ratio = g/total_len
t_ratio = t/total_len

# Output
print(f'Count: {len(sizes)}')
print(f'Total: {total_len}')
print(f'Min: {min}')
print(f'Max: {max}')
print(f'Mean: {total_len/len(sizes):.1f}')
print(f'NTs: {a_ratio:.3f} {c_ratio:.3f} {g_ratio:.3f} {t_ratio:.3f}')
#print(sizes[int(len(sizes)/2)])

"""
python3 fasta_stats.py transcripts.fasta.gz
Count: 232
Total: 278793
Min: 603
Max: 1991
Mean: 1201.7
NTs: 0.291 0.218 0.210 0.281
"""
