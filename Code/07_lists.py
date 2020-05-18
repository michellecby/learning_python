#!/usr/bin/env python3

import math
import random

# Move the triple quotes downward to uncover each segment of code


# Sometimes we want to work with multiple variables at the same time
# A convenient way to do this is a 'tuple'
# Any time you have comma separated values in parentheses, you have a 'tuple'
# Tuples can mix data types

tup = (1, 2.0, 'three')

# You can access any element with square bracket context
# Counting starts at zero, not 1

print(tup[0])
print(tup[2])

# Tuples are 'immutable', meaning you can't change their contents

#tup[2] = 3 # uncomment this line to try it out

# Lists are like tuples, but are mutable
# Create lists with square brackets instead of parentheses

list = [1, 2.0, 'three']
print(list)

# You can change the content of each item in a list

list[2] = 3
print(list)

# You can get a slice of a list, which should look familiar

print(list[0:2])

# You can change the length of a list with append(), pop(), and insert()

list.append('four') # adds 'four' to the end of a list
print(list)

last = list.pop() # removes the last element of a list and returns it
print(last, list)

list.insert(2, 'ok')  # inserts 'ok' at position 2
print(list)

# Many lists contain numbers, like the following probability distribution

p = [0.2, 0.1, 0.4, 0.3]

# Note the two ways of iterating over a list

sum1 = 0
for x in p:
	sum1 += x

sum2 = 0
for i in range(len(p)):
	sum2 += p[i]

print(sum1, sum2)

# If all the elements of a list are numbers, it's easy to sort

numbers = [5, 1, 15, 2]
numbers.sort()
print(numbers)

# If all the elements are strings, it's also easy to sort

words = ['yo', 'hello', 'world', 'arrg', 'a', 'a2', 'A2', '2a']
words.sort()
print(words)

# To convert from strings to lists, use split()

text = 'default split      is     any     spaces'
cols = text.split()
print(cols)

line = '1,2,3'
csv = line.split(',') # comma separated values
print(csv) # note the quotes, the numeric values are actually strings


# If you want to split every letter of a string, use list()

# To convert from lists to strings, use join()
letters = list('ACGT')
print(letters)

seq = ['G', 'A', 'A', 'T', 'T', 'C']
s = ''.join(seq)
print(s, seq)

# The command line is a list called sys.argv
# Try this:
# 	python3 07_lists.py a b c 1 2 3

import sys
print(sys.argv)
