#!/usr/bin/env python3

import math

# Move the triple quotes downward to uncover each segment of code


# We have been using functions for a while, such as print(list)
# Sometimes the variable preceeds the function, such as list.append()
# In this section, we will learn how to write new functions, like print()

# A function is created with the `def` keyword followed by the name
# Functions are usually declared at the top of the file
# You can't call a function until it has been 'seen' by python
# Functions usually have parameters, but this next one does not

def greeting():
	print('hey')

for i in range(3): greeting()

# Function parameters go in parentheses
# In the next function, there is one parameter: n

def stirling(n):
	if n == 0: return 0
	return 0.5 * math.log(math.tau) + (n + 0.5) * math.log(n) \
        - n + 1/(12 * n) - 1 / (360 * (n**3))

for i in range(6):
	print(i, 'factorial is approximately', math.e ** stirling(i))

# Functions can have several parameters, as a tuple

def poisson(n, k):
	p = k ** n * math.e ** -k / math.e ** stirling(n)
	return p

for n in range(10):
	print('poisson prob of %d with expectation 4 %.4f' % (n, poisson(n, 4)))

# Functions can return several return values, as a tuple

def quadratic(a, b, c):
	d = b ** 2 - 4 * a * c
	assert(d > 0)
	x1 = (-b + math.sqrt(d)) / (2 * a)
	x2 = (-b - math.sqrt(d)) / (2 * a)
	return x1, x2

print(quadratic(2, 5, 3))

# Functions can take lists as parameters

def max(values):
	assert(len(values) > 0)
	max = values[0]
	for v in values:
		if v > max: max = v
	return max

data = [-1, 10, 2, 0, 3.14]
print(max(data))

# Functions can return lists too

def min_max(values):
	assert(len(values) > 0)
	min = values[0]
	max = values[0]
	for v in values:
		if v > max: max = v
		if v < min: min = v
	return [min, max]

print(min_max(data))



# Here's a really useful function for reading fasta files
# Some parts of the following code may be unfamiliar, don't stress about it!
# Just copy-paste read_fasta into the top of your program
# Files may be compressed (.gz suffix) or piped in (use - for filename)

import gzip

def read_fasta(filename):
	name = None
	seqs = []
	
	fp = None
	if filename == '-':
		fp = sys.stdin
	elif filename.endswith('.gz'):
		fp = gzip.open(filename, 'rt')
	else:
		fp = open(filename)

	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield(name, seq)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield(name, ''.join(seqs))
	fp.close()



for name, seq in read_fasta('proteins.fasta.gz'):
	print(name, len(seq))

# Note that copy-pasting functions between programs is a deplorable practice
# We will soon learn better methods of code re-use


