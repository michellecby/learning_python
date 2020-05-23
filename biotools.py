#!/usr/bin/env python3

import sys
import gzip
import random
import math

# Index
# 1. read_fasta: read fasta file
# 2. gc: count gc ratio
# 3. randseq: generate random sequence with specific gc ratio
# 4. hydro: calculate hydropobicity score in kd, is, or os (dict)
# 5. kdscore: calculate kd score in a window (elif)
# 6. translation: translate nt to aa (dict)
# 7. trans_elif: translate nt to aa (elif)
# 8. skew: calculate gc-skew(g and c are over or under abundant)
# 9. entropy: calculate Shannon's entropy in the sequence (dict)

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

kdscale = { 'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,
	   'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
	   'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
	   'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }

isscale = { 'I':-0.31,'L':-0.56,'F':-1.13,'V':0.07,'M':-0.23,
	   'P':0.45,'W':-1.85,'H': 0.17,'T':0.14,'E':-0.01,
	   'Q':0.58,'C':-0.24,'Y':-0.94,'A':0.17,'S': 0.13,
	   'N':0.42,'D':-0.07,'R': 0.81,'G':0.01,'K': 0.99 }

osscale = { 'I':-1.12,'L':-1.25,'F':-1.71,'V':-0.46,'M':-0.67,
	   'P':0.14,'W':-2.09,'H': 0.11,'T':0.25,'E': 0.11,
	   'Q':0.77,'C':-0.02,'Y':-0.71,'A':0.50,'S':0.46,
	   'N':0.85,'D': 0.43,'R': 1.81,'G':1.15,'K':2.80 }

def hydro(protein, method):
	assert(method == 'kd' or method == 'is' or method == 'os')
	hydro_score = 0
	if method == 'kd':
		scale = kdscale
	elif method == 'is':
		scale = isscale
	elif method == 'os':
		scale = osscale

	for aa in protein:
		if aa == '*': continue
		hydro_score += scale[aa]
	return hydro_score

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

trans = {'ATG':'M', 'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':"L",\
 'CTC':'L', 'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'GTT':'V',\
 'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',\
 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T', 'ACC':'T', 'ACA':'T',\
 'ACG':'T', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'TAT':'Y', 'TAC':'Y',\
 'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'AAT':'N', 'AAC':'N', 'AAA':'K',\
 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'TGT':'C', 'TGC':'C',\
 'TGG':'W', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGT':'S', 'AGC':'S',\
 'AGA':'R', 'AGG':'R', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 'TGA':'*',\
 'TAA':'*', 'TAG':'*'}
def translation(sequence):
	prot = []
	assert(len(sequence) % 3 == 0)
	for i in range(0, len(sequence)-3+1, 3):
		codon = sequence[i:i+3]
		if codon not in trans: prot.append('X')
		prot.append(trans[codon])
	return ''.join(prot)

def trans_elif(sequence):
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
	count = {'A':0, 'C':0, 'G':0, 'T':0}
	for nt in win:
		if   nt == 'A': count['A'] += 1
		elif nt == 'C': count['C'] += 1
		elif nt == 'G': count['G'] += 1
		elif nt == 'T': count['T'] += 1

	total = count['A'] + count['C'] + count['G'] + count['T']
	h = 0
	pa, pc, pg, pt = count['A']/total, count['C']/total, count['G']/total, count['T']/total

	if count['A'] != 0: h -= pa * math.log2(pa)
	if count['C'] != 0: h -= pc * math.log2(pc)
	if count['G'] != 0: h -= pg * math.log2(pg)
	if count['T'] != 0: h -= pt * math.log2(pt)
	return h

"""
def entropy(seq):
	h = 0
	a, c, g, t = 0, 0, 0, 0
	for nt in seq:
		if   nt == 'A': a += 1
		elif nt == 'C': c += 1
		elif nt == 'G': g += 1
		elif nt == 'T': t += 1
	total = a + c + g + t
	pa, pc, pg, pt = a/total, c/total, g/total, t/total
	if count == 0:
		return None
	if a > 0: h -= pa * math.log2(pa)
	if c > 0: h -= pc * math.log2(pc)
	if g > 0: h -= pg * math.log2(pg)
	if t > 0: h -= pt * math.log2(pt)
	return h
"""
