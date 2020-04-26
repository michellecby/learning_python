#!/usr/bin/env python3

# Write a program that computes the GC fraction of a DNA sequence in a window
# Window size is 11 nt
# Step size is 5 nt
# Output with 4 significant figures using whichever method you prefer
# Use no nested loops
# Describe the pros/cons of this algorith vs. nested loops

seq = 'ACGACGCAGGAGGAGAGTTTCAGAGATCACGAATACATCCATATTACCCAGAGAGAG'
w = 11
s = 5

gc = 0
for i in range(0, w):
    if seq[i] == "G" or seq[i] == "C":
        gc += 1
print(f'0 {seq[0:w]} {gc/w:.4f}')

for i in range(s, len(seq)-w+1, s):
    prev = seq[i-s:i]
    next = seq[i+w-s:i+w]
    for p in prev:
        if p == "G" or p == "C":
            gc -= 1
    for n in next:
        if n == "G" or n == "C":
            gc += 1
    print(f'{i} {seq[i:i+w]} {gc/w:.4f}')

"""
0 ACGACGCAGGA 0.6364
5 GCAGGAGGAGA 0.6364
10 AGGAGAGTTTC 0.4545
15 AGTTTCAGAGA 0.3636
20 CAGAGATCACG 0.5455
25 ATCACGAATAC 0.3636
30 GAATACATCCA 0.3636
35 CATCCATATTA 0.2727
40 ATATTACCCAG 0.3636
45 ACCCAGAGAGA 0.5455
"""
