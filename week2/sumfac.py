#!/usr/bin/env python3

# Write a program that computes the running sum from 1..n
# Also, have it compute the factorial of n while you're at it
# No, you may not use math.factorial()
# Use the same loop for both calculations

n = 5
num_sum=0
num_fac=1
for i in range(1, n+1):
    num_sum+=i
    num_fac=num_fac*i

print(n, num_sum, num_fac)
"""
5 15 120
"""
