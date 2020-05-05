#!/usr/bin/env python3

# Print out all the codons for the sequence below in reading frame 1
# Use a 'for' loop

dna = 'ATAGCGAATATCTCTCATGAGAGGGAA'

for i in range(0, len(dna), 3):
    #print(dna[i]+dna[i+1]+dna[i+2])
    codon=(dna[i:i+3])
    print(codon)

"""
for i in range(3):
    print('frame',f)
    for i in range(0, len(dna)-k+1, 3):
        kmer=dna[i:i+k]
        print(kmer)
"""

"""
ATA
GCG
AAT
ATC
TCT
CAT
GAG
AGG
GAA
"""
