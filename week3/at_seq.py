#!/usr/bin/env python3

import random
#random.seed(1) # comment-out this line to change sequence each time

length = 30
dna = ''
w_count =0

for i in range(length):
    r = random.random() #generate a random number from 0 <= r < 1
    if r <= 0.6:
        w = random.randint(0,1)
        w_count += 1
        if w == 0:
            dna += "A"
        else:
            dna += "T"
    else:
        s = random.randint(0,1)
        if s == 0:
            dna += "G"
        else:
            dna += "C"
print(length, w_count/length, dna)


# 50% AT ratio
"""
alpha = "AGCT"
for i in range(length):
    nt = random.choice(alpha)
    dna += nt

for j in dna:
    if j == "G" or j == "C":
        W_count += 1
print(length, W_count/length, dna)
"""


# Write a program that stores random DNA sequence in a string
# The sequence should be 30 nt long
# On average, the sequence should be 60% AT
# Calculate the actual AT fraction while generating the sequence
# Report the length, AT fraction, and sequence



"""
30 0.6666666666666666 ATTACCGTAATCTACTATTAAGTCACAACC
"""
