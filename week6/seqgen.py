import argparse
import random
import biotools
# setup
parser = argparse.ArgumentParser(
	description='Generate random FASTA file.')
# required arguments
parser.add_argument('--count', required=True, type=int,
	metavar='<int>', help='number of sequences')
parser.add_argument('--min', required=True, type=int,
	metavar='<int>', help='minimum sequence length')
parser.add_argument('--max', required=True, type=int,
	metavar='<int>', help='maximum sequence length')

parser.add_argument('--gc', required=True, type=float,
	metavar='<float>', help='the gc composition')

# optional arguments with default parameters
parser.add_argument('--prefix', required=False, type=str, default='seq',
	metavar='<str>', help='optional string argument [%(default)s]')

# switches
parser.add_argument('--verbose', action='store_true',
	help='on/off switch')
# finalization
arg = parser.parse_args()

for i in range(arg.count):
    print(f'>{arg.prefix}_{i}')
    seq_len = random.randint(arg.min, arg.max)
    print(biotools.randseq(seq_len, arg.gc))
