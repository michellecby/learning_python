#!/usr/bin/env python3

# Write a program that computes hydrophobicity in a window
# Let the user choose the method (see below)
# https://en.wikipedia.org/wiki/Hydrophilicity_plot
# https://en.wikipedia.org/wiki/Hydrophobicity_scales
import argparse
import biotools

parser = argparse.ArgumentParser(
	description='Hydrophobicity calculator')
# required arguments
parser.add_argument('--input', required=True, type=str,
	metavar='<path>', help='file input')
parser.add_argument('--window', required=False, type=int, default = 15,
	metavar='<int>', help='Window size, default [%(default)i]')
parser.add_argument('--method', required=False, type=str, default = 'kd',
	metavar='<str>', help='Choose one calculation method: kd (Kyte-Doolittle scale), is (Interface scale), os (Octnal scale). Default: [%(default)s]')

# finalization
arg = parser.parse_args()

assert(arg.method == 'kd' or arg.method == 'is' or arg.method == 'os')
for name, seq in biotools.read_fasta(arg.input):
	with open(name+'_'+arg.method+'.tsv', "w") as fp:
		for i in range(0, len(seq)-arg.window+1, arg.window):
			peptide = seq[i:i+arg.window]
			fp.write(f'{i}\t{biotools.hydro(peptide, arg.method)/len(peptide):.4f}\n')


"""
python3 hydrophobicity.py --input proteins.fasta.gz --window 11 --method kd
"""
