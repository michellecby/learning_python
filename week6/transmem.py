#!/usr/bin/env python3

import biotools
import argparse

# Use argparse
# Write a program that translates an mRNA
# Assume the protein encoded is the longest ORF

# setup
parser = argparse.ArgumentParser(
	description='Predict transmembrane proteins.')
# required arguments
parser.add_argument('--file', required=True, type=str,
	metavar='<path>', help='Protein file')
# optional arguments with default parameters
parser.add_argument('--wind1', required=False, type=int, default=8,
	metavar='<int>', help='length of signal peptide region [%(default)i]')
parser.add_argument('--wind2', required=False, type=int, default=11,
	metavar='<int>', help='length of hydrophobic region [%(default)i]')
parser.add_argument('--kd1', required=False, type=float, default=2.5,
	metavar='<float>', help='kd value for signal peptide region [%(default)f]')
parser.add_argument('--kd2', required=False, type=float, default=2.0,
	metavar='<float>', help='kd value for hydrophobic region [%(default)f]')

# finalization
arg = parser.parse_args()


def hydro(seq, wind, threhd):
	for i in range(len(seq) - wind +1):
		peptide = seq[i:i + wind]
		if biotools.kdscore(peptide) > threhd and 'P' not in peptide:
			return True
	return False

for name, seq in biotools.read_fasta(arg.file):
	nterm = seq[:30]
	cterm = seq[30:]
	if hydro(nterm, arg.wind1, arg.kd1) and hydro(cterm, arg.wind2, arg.kd2):
		print(name)

"""
python3 transmem.py --file ../week5/proteins.fasta.gz

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
