#!/usr/bin/env python3

# Write a program that computes the GC% of a DNA sequence
# Format the output for 2 decimal places
# Use all three formatting methods

dna = 'ACAGAGCCAGCAGATATACAGCAGATACTAT' # feel free to change
s_count = 0
for nt in dna:
    if nt == "G" or nt == "C":
        s_count += 1

s_frac = s_count / len(dna)

# Method 1: %
print('%.2f' %(s_frac))

# Method 2:{}.format
print('{:.2f}'.format(s_frac))

#Method 3:f{}
print(f'{s_frac:.2f}')


"""
0.42
0.42
0.42
"""
