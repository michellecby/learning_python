#!/usr/bin/env python3

# Write a program that prints the reverse-complement of a DNA sequence
# You must use a loop and conditional

dna = 'ACTGAAAAAAAAAAA'
rev_comp = ''
for nt in dna:
    if nt == "A":
        rev_comp = "T" + rev_comp
    elif nt == "T":
        rev_comp = "A" + rev_comp
    elif nt == "C":
        rev_comp = "G" + rev_comp
    else:
        rev_comp = "C" + rev_comp
print(dna)
print(rev_comp)

"""
TTTTTTTTTTTCAGT
"""
