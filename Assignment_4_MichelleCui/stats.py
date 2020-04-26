#!/usr/bin/env python3

from math import sqrt
import fileinput

# Write a program that computes typical stats
# Count, Min, Max, Mean, Std. Dev, Median
# No, you cannot import any other modules!

lst = []
for line in fileinput.input():
    if line.startswith('#'): continue
    num = float(line.rstrip())
    lst.append(num)
lst.sort()
#print(lst)

count = len(lst)
min = lst[0]
max = lst[-1]
sum = 0
for i in lst:
    sum += i
mean = sum/count

var = 0
for j in lst:
    var += (j - mean) * (j - mean)
std = sqrt(var/count)

m = int(count/2)-1
if count%2 != 0:
    median = lst[round(m)+1]
else:
    median = (lst[m]+lst[m+1])/2

# output
print(f'Count: {count}')
print(f'Minimum: {min}')
print(f'Maximum: {max}')
print(f'Mean: {mean}')
print(f'Std. dev: {std:.5}')
print(f'Meidan: {median}')


"""
python3 stats.py numbers.txt
Count: 10
Minimum: -1.0
Maximum: 256.0
Mean: 29.147789999999997
Std. dev: 75.777
Median 2.35914
"""
