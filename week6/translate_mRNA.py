#!/usr/bin/env python3

import biotools
import argparse

# Use argparse
# Write a program that translates an mRNA
# Assume the protein encoded is the longest ORF
# setup
parser = argparse.ArgumentParser(
	description='Reports the longest coding sequence.')
# required arguments
parser.add_argument('--file', required=True, type=str,
	metavar='<path>', help='Sequence file.')

# finalization
arg = parser.parse_args()

def orfseq(sequence):
	assert(len(sequence) > 0)

	cds = None
	max_len = 0
	atgs=[]
	for i in range(len(sequence)-2):
		codon = sequence[i:i+3]
		if codon == "ATG":
			atgs.append(i)
	for atg in atgs:
		for i in range(atg, len(sequence), 3):
			codon = sequence[i:i+3]
			if codon == "TAA" or codon == "TGA" or codon == "TAG":
				stop = i
				orf = stop - atg + 1
				if orf > max_len:
					max_len = orf
					cds = sequence[atg:stop+3]
				break
		# if stop != None:
			#orf = stop - atg + 3
			#if orf > max_len:
				#max_len = orf
				#cds = sequence[atg:atg + orf]
	return cds

for name, seq in biotools.read_fasta(arg.file):
	print(f'>{name}')
	print(biotools.translation(orfseq(seq)))

"""
python3 translate_mRNA.py --file ../week5/transcripts.fasta.gz
>CBG00001.1
MTFCENKNLPKPPSDRCQVVVISILSMILDFYLKYNPDKHWAHLFYGASPILEILVIFGMLANSVYGNKLAMFACVLDLVSGVFCLLTLPVISVAENATGVRLHLPYISTFHSQFSFQVSTPVDLFYVATFLGFVSTILILLFLILDALKFMKLRKLRNEDLEKEKKMNPIEKV*
>CBG00006.1
MNGVEKVNKYFDIKDKRDFLYHFGFGVDTLDIKAVFGDTKFVCTGGSPGRFKLYAEWFAKETSIPCSENLSRSDRFVIYKTGPVCWINHGMGTPSLSIMLVESFKLMHHAGVKNPTFIRLGTSGGVGVPPGTVVVSTGAMNAELGDTYVQVIAGKRIERPTQLDATLREALCAVGKEKNIPVETGKTMCADDFYEGQMRLDGYFCDYEEEDKYAFLRKLNSLGVRNIEMESTCFASFTCRAGFPSAIVCVTLLNRMDGDQVQIDKEKYIEYEERPFRLVTAYIRQQTGV*
etc.
"""
