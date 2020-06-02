#!/usr/bin/env python3

import argparse
import biotools

parser = argparse.ArgumentParser(
	description='Position weight matrix finder.')
parser.add_argument('--dna', required=True, type=str,
	metavar='<str>', help='DNA FASTA file')
parser.add_argument('--pwm', required=True, type=str,
	metavar='<str>', help='Position weight matrix transfac file')
parser.add_argument('--threshold', required=True, type=float,
	metavar='<float>', help='probability threshold for finding motif')
arg = parser.parse_args()

def read_transfac(filename):
	pwm = []
	assert(filename.endswith('.transfac'))

	fp = None
	if filename == '-':
		fp = sys.stdin
	else:
		fp = open(filename)

	for line in fp.readlines()[6: ]:
		line = line.rstrip()
		line = line.split('\t')
		if len(line) > 0:
			if line[0] != 'XX':
				pwm.append({'A': line[1], 'C': line[2], 'G': line[3], 'T': line[4]})
			else:
				break
	return pwm



def pwm_prob (sequence, thre):

	total = 0
	for key, value in pos_weight[0].items():
		total += float(value)

	win = len(pos_weight)
	pos = {}
	for i in range(len(sequence)-win+1):
		kmer = sequence[i:i+win]
		prob = 1
		for j in range(len(kmer)):
			prob *= float(pos_weight[j][str(kmer[j])])/total
		if prob > thre:
			pos[i]=prob
	return pos


pos_weight = read_transfac(arg.pwm)

for name, seq in biotools.read_fasta(arg.dna):
	name = name.split(' ')
	for position, probability in pwm_prob(seq, arg.threshold).items():
		print(f'{name[0]}\t{position}\t{seq[position:position+len(pos_weight)]}\t{probability:.4f}')


"""
python3 pwm_search.py --dna sars-cov-2.fa.gz --pwm MA0036.1.transfac --threshold 0.01
"""
