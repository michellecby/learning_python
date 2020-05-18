#!/usr/bin/env python3

import gzip
import sys

# Write a program that reports the longest coding sequence
# Translate the sequence to amino acids
# Use a generator function somewhere in your program
# Check both strands
# Command line:
#	python3 longest_cds.py transcripts.fasta.gz


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

def translation(sequence):
	aa = []
	for i in range(0, len(sequence)-3+1, 3):
		codon = sequence[i:i+3]
		if codon == 'ATG':
			aa.append('M')
		elif codon == 'TTT' or codon == 'TTC':
			aa.append('F')
		elif codon == 'TTA' or codon == 'TTG':
			aa.append('L')
		elif codon == 'CTT' or codon == 'CTC' or codon == 'CTA' or codon == 'CTG':
			aa.append('L')
		elif codon == 'ATT' or codon == 'ATC' or codon == 'ATA':
			aa.append('I')
		elif codon == 'GTT' or codon == 'GTC' or codon == 'GTA' or codon == 'GTG':
			aa.append('V')
		elif codon == 'TCT' or codon == 'TCC' or codon == 'TCA' or codon == 'TCG':
			aa.append('S')
		elif codon == 'CCT' or codon == 'CCC' or codon == 'CCA' or codon == 'CCG':
			aa.append('P')
		elif codon == 'ACT' or codon == 'ACC' or codon == 'ACA' or codon == 'ACG':
			aa.append('T')
		elif codon == 'GCT' or codon == 'GCC' or codon == 'GCA' or codon == 'GCG':
			aa.append('A')
		elif codon == 'TAT' or codon == 'TAC':
			aa.append('Y')
		elif codon == 'CAT' or codon == 'CAC':
			aa.append('H')
		elif codon == 'CAA' or codon == 'CAG':
			aa.append('Q')
		elif codon == 'AAT' or codon == 'AAC':
			aa.append('N')
		elif codon == 'AAA' or codon == 'AAG':
			aa.append('K')
		elif codon == 'GAT' or codon == 'GAC':
			aa.append('D')
		elif codon == 'GAA' or codon == 'GAG':
			aa.append('E')
		elif codon == 'TGT' or codon == 'TGC':
			aa.append('C')
		elif codon == 'TGG':
			aa.append('W')
		elif codon == 'CGT' or codon == 'CGC' or codon == 'CGA' or codon == 'CGG':
			aa.append('R')
		elif codon == 'AGT' or codon == 'AGC':
			aa.append('S')
		elif codon == 'AGA' or codon == 'AGG':
			aa.append('R')
		elif codon == 'GGT' or codon == 'GGC' or codon == 'GGA' or codon == 'GGG':
			aa.append('G')
	return ''.join(aa)


protein_file = sys.argv[1]

def orfseq(sequence):
	cds = None
	max_len = 0
	# find all ATG
	atgs=[]
	for i in range(len(sequence)-2):
		codon = sequence[i:i+3]
		if codon == "ATG":
			atgs.append(i)
	# for every ATG
	for atg in atgs:
		#find closest stop
		#stop = None
		for i in range(atg, len(sequence), 3):
			codon = sequence[i:i+3]
			if codon == "TAA" or codon == "TGA" or codon == "TAG":
				stop = i
				orf = stop - atg + 1
				if orf > max_len:
					max_len = orf
					cds = sequence[atg:stop+3]
				break
	return cds
		# record the len
		#if stop != None:

for name, seq in read_fasta(protein_file):
	print(f'>{name}')
	print(translation(orfseq(seq)))

"""
>CBG00001.1
MTFCENKNLPKPPSDRCQVVVISILSMILDFYLKYNPDKHWAHLFYGASPILEILVIFGMLANSVYGNKLAMFACVLDLVSGVFCLLTLPVISVAENATGVRLHLPYISTFHSQFSFQVSTPVDLFYVATFLGFVSTILILLFLILDALKFMKLRKLRNEDLEKEKKMNPIEKV*
>CBG00006.1
MNGVEKVNKYFDIKDKRDFLYHFGFGVDTLDIKAVFGDTKFVCTGGSPGRFKLYAEWFAKETSIPCSENLSRSDRFVIYKTGPVCWINHGMGTPSLSIMLVESFKLMHHAGVKNPTFIRLGTSGGVGVPPGTVVVSTGAMNAELGDTYVQVIAGKRIERPTQLDATLREALCAVGKEKNIPVETGKTMCADDFYEGQMRLDGYFCDYEEEDKYAFLRKLNSLGVRNIEMESTCFASFTCRAGFPSAIVCVTLLNRMDGDQVQIDKEKYIEYEERPFRLVTAYIRQQTGV*
etc.
"""
