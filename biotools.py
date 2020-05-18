#!/usr/bin/env python3

import sys
import gzip
import random
import math

# Index
# 1. read_fasta: read fasta file
# 2. gc: count gc ratio
# 3. randseq: generate random sequence with specific gc ratio
# 4. kdscore: calculate kd score in a window
# 5. translation: translate nucleotide sequence to amino aci sequence
# 6. skew: calculate gc-skew(g and c are over or under abundant)
# 7. entropy: calculate Shannon's entropy in the sequence

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

def gc(seq):
	count = 0
	for nt in seq:
		if nt == 'G' or nt == 'C':
			count += 1
	return count / len(seq)

def randseq(l, gc):
	seq = []
	for i in range(l):
		r = random.random()
		if r < gc:
			r = random.random()
			if r > 0.5: seq.append("G")
			else:       seq.append("C")
		else:
			r = random.random()
			if r > 0.5: seq.append("A")
			else:       seq.append("T")
	return ''.join(seq)

def kdscore(window):
	kd = 0
	for aa in window:
		if aa == 'I':   kd += 4.5
		elif aa == 'V': kd += 4.2
		elif aa == 'L': kd += 3.8
		elif aa == 'F': kd += 2.8
		elif aa == 'C': kd += 2.5
		elif aa == 'M': kd += 1.9
		elif aa == 'A': kd += 1.8
		elif aa == 'G': kd -= 0.4
		elif aa == 'T': kd -= 0.7
		elif aa == 'S': kd -= 0.8
		elif aa == 'W': kd -= 0.9
		elif aa == 'Y': kd -= 1.3
		elif aa == 'P': kd -= 1.6
		elif aa == 'H': kd -= 3.2
		elif aa == 'E': kd -= 3.5
		elif aa == 'Q': kd -= 3.5
		elif aa == 'D': kd -= 3.5
		elif aa == 'N': kd -= 3.5
		elif aa == 'K': kd -= 3.9
		elif aa == 'R': kd -= 4.5
	#print(window, kd)
	return kd/len(window)

def translation(sequence):
	prot = []
	assert(len(sequence) % 3 == 0)
	for i in range(0, len(sequence)-3+1, 3):
		codon = sequence[i:i+3]
		if codon == 'ATG':
			prot.append('M')
		elif codon == 'TTT' or codon == 'TTC':
			prot.append('F')
		elif codon == 'TTA' or codon == 'TTG':
			prot.append('L')
		elif codon == 'CTT' or codon == 'CTC' or codon == 'CTA' or codon == 'CTG':
			prot.append('L')
		elif codon == 'ATT' or codon == 'ATC' or codon == 'ATA':
			prot.append('I')
		elif codon == 'GTT' or codon == 'GTC' or codon == 'GTA' or codon == 'GTG':
			prot.append('V')
		elif codon == 'TCT' or codon == 'TCC' or codon == 'TCA' or codon == 'TCG':
			prot.append('S')
		elif codon == 'CCT' or codon == 'CCC' or codon == 'CCA' or codon == 'CCG':
			prot.append('P')
		elif codon == 'ACT' or codon == 'ACC' or codon == 'ACA' or codon == 'ACG':
			prot.append('T')
		elif codon == 'GCT' or codon == 'GCC' or codon == 'GCA' or codon == 'GCG':
			prot.append('A')
		elif codon == 'TAT' or codon == 'TAC':
			prot.append('Y')
		elif codon == 'CAT' or codon == 'CAC':
			prot.append('H')
		elif codon == 'CAA' or codon == 'CAG':
			prot.append('Q')
		elif codon == 'AAT' or codon == 'AAC':
			prot.append('N')
		elif codon == 'AAA' or codon == 'AAG':
			prot.append('K')
		elif codon == 'GAT' or codon == 'GAC':
			prot.append('D')
		elif codon == 'GAA' or codon == 'GAG':
			prot.append('E')
		elif codon == 'TGT' or codon == 'TGC':
			prot.append('C')
		elif codon == 'TGG':
			prot.append('W')
		elif codon == 'CGT' or codon == 'CGC' or codon == 'CGA' or codon == 'CGG':
			prot.append('R')
		elif codon == 'AGT' or codon == 'AGC':
			prot.append('S')
		elif codon == 'AGA' or codon == 'AGG':
			prot.append('R')
		elif codon == 'GGT' or codon == 'GGC' or codon == 'GGA' or codon == 'GGG':
			prot.append('G')
		elif codon == 'TGA' or codon == 'TAA' or codon == 'TAG':
			prot.append('*')
		else: prot.append('X')
	return ''.join(prot)

def skew(seq):
	#(g-c)/(g+c)
	g = 0
	c = 0
	for nt in seq:
		if nt == 'G':
			g += 1
		elif nt == 'C':
			c += 1
	return (g - c)/(g + c)


def entropy(seq):
	h = 0
	a = 0
	c = 0
	g = 0
	t = 0
	for nt in seq:
		if   nt == 'A': a += 1
		elif nt == 'C': c += 1
		elif nt == 'G': g += 1
		elif nt == 'T': t += 1
	count = a + c + g + t
	if count == 0:
		return None
	if a > 0: h -= (a/count) * math.log2(a/count)
	if c > 0: h -= (c/count) * math.log2(c/count)
	if g > 0: h -= (g/count) * math.log2(g/count)
	if t > 0: h -= (t/count) * math.log2(t/count)
	return h
