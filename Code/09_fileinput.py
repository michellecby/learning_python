#!/usr/bin/env python3

import fileinput
import math
import sys

# Move the triple quotes downward to uncover each segment of code



# We have been using hard-coded data in our programs like this

dna = 'CATGGGCACGCATTACAGT'

# This is a terrible practice, and we won't do it anymore
# Data should come from outside the program
# Fortunately, python makes that very simple with the fileinput module
# It reads every line of every file on the command line

# Run this:
# 	python3 09_fileinput.py numbers.txt

data = []
for line in fileinput.input():
	#if line[0] == '#': continue # skip over comments
	if line.startswith('#'): continue # same as above
	line = line.rstrip() # remove newline (return character), often useful
	data.append(float(line)) # store the data


# Now all of your data is in a list and you can do stuff with it
data.sort()
print('range = {} to {}'.format(data[0], data[-1]))
