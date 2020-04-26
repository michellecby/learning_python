#!/usr/bin/env python3

import fileinput

# Write a program that computes typical sequence stats
# No, you cannot import any other modules!
# Use rand_seq to generate the sequences
# Expected output is shown below

seq_length = []
letter_count = 0
a_count = 0
c_count = 0
g_count = 0
t_count = 0
for line in fileinput.input():
    if line.startswith('>'): continue
    seq = line.rstrip()
    letter_count += len(seq)
    seq_length.append(len(seq))
    for nt in seq:
        if nt == "A":
            a_count += 1
        elif nt == "C":
            c_count += 1
        elif nt == "G":
            g_count += 1
        else:
            t_count += 1

seq_length.sort()
min = seq_length[0]
max = seq_length[1]

threshold = letter_count/2
sum = 0
n50 = None
for x in seq_length:
    sum += x
    n50 = x
    if sum >= threshold:
        break



#output
print(f'Number of sequences: {len(seq_length)}')
print(f'Number of letters: {letter_count}')
print(f'Minimum length: {min}')
print(f'Maximum lenth: {max}')
print(f'N50: {n50}')
print(f'Composition: A={a_count/letter_count:.3} C={c_count/letter_count:.3} G={g_count/letter_count:.3} T={t_count/letter_count:.3}')

"""
python3 rand_seq.py 100 100 100000 0.35 | python3 seqstats.py
Number of sequences: 100
Number of letters: 4957689
Minimum length: 219
Maximum length: 99853
N50: 67081
Composition: A=0.325 C=0.175 G=0.175 T=0.325
"""
