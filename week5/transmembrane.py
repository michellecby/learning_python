#!/usr/bin/env python3

import gzip
import sys

# Write a program that predicts if a protein is trans-membrane
# Trans-membrane proteins have the following properties
#	Signal peptide: https://en.wikipedia.org/wiki/Signal_peptide
#	Hydrophobic regions(s): https://en.wikipedia.org/wiki/Transmembrane_protein
#	No prolines (alpha helix) in hydrophobicity part
# Hydrophobicity is measued via Kyte-Dolittle
#	https://en.wikipedia.org/wiki/Hydrophilicity_plot
# For our purposes:
#	Signal peptide is 8 aa long, KD > 2.5, first 30 aa
#	Hydrophobic region is 11 aa long, KD > 2.0, after 30 aa
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

def hydro(seq, wind, threhd):
	for i in range(len(seq) - wind +1):
		peptide = seq[i:i + wind]
		if kdscore(peptide) > threhd and 'P' not in peptide:
			return True
	return False

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

for name, seq in read_fasta(protein_file):
	nterm = seq[:30]
	cterm = seq[30:]
	if hydro(nterm, 8, 2.5) and hydro(cterm, 11, 2.0):
		#if 'P' not in nterm and 'P' not in cterm:
		print(name)



"""
python3 transmembrane.py proteins.fasta.gz

18w
Dtg
Krn
Lac
Mcr
PRY
Pxt
Pzl
QC
Ror
S1P
S2P
Spt
apn
bai
bdl
bou
bug
cue
drd
ft
grk
knk
ksh
m
nac
ort
rk
smo
thw
tsg
waw
zye
"""
