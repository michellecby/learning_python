#!/usr/bin/env python3

import argparse
import biotools

# Write a program that predicts if a protein is trans-membrane
# Trans-membrane proteins have the following properties
#	Signal peptide: https://en.wikipedia.org/wiki/Signal_peptide
#	Hydrophobic regions(s): https://en.wikipedia.org/wiki/Transmembrane_protein
#	No prolines (alpha helix)
# Hydrophobicity is measued via Kyte-Dolittle
#	https://en.wikipedia.org/wiki/Hydrophilicity_plot
# For our purposes:
#	Signal peptide is 8 aa long, KD > 2.5, first 30 aa
#	Hydrophobic region is 11 aa long, KD > 2.0, after 30 aa

# Use a dictionary for the K-D values, and store in biotools
# Use argparse for command line
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


def hydro_check(seq, wind, threhd):
	for i in range(len(seq) - wind +1):
		peptide = seq[i:i + wind]
		if (biotools.hydro(peptide, 'kd')/len(peptide)) > threhd and 'P' not in peptide:
			return True
	return False

for name, seq in biotools.read_fasta(arg.file):
	nterm = seq[:30]
	cterm = seq[30:]
	if hydro_check(nterm, arg.wind1, arg.kd1) and hydro_check(cterm, arg.wind2, arg.kd2):
		print(name)
"""
python3 transmembrane.py --file proteins.fasta.gz
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
