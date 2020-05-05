#!/usr/bin/env python3

# Move the triple quotes downward to uncover each segment of code


# The most common way to get data into a program is by reading a file
# Before reading a file you must open it first
# open() associates a file with a 'file pointer' variable
# A file pointer is its own kind of variable type

print('\n--- type ---')
with open('numbers.txt', 'r') as fp:
	print(type(fp))

# In the code above, the 'with' statment creates block structure
# The file is automatically closed after leaving the block
# The 'r' means 'read' and is the default, so it can be omitted (as below)

# You can read from the file pointer several ways
# The read() method returns the entire contents of a file

print('\n--- read() ---')
with open('numbers.txt') as fp:
	stuff = fp.read()
print(stuff)

# The readline() method returns a single line

print('\n--- readline() ---')
with open('numbers.txt') as fp:
	line = fp.readline()
print(line)

# The readlines() method returns all of the lines as an array of strings
# That means you can iterate over the contents with a 'for' loop

print('\n--- readlines() ---')
numbers = []
with open('numbers.txt') as fp:
	for line in fp.readlines():
		if line[0] != '#':
			numbers.append(float(line))
print(numbers)

# Biological data files can be immense, and are often compressed with gzip
# You can read from compressed files directly without uncompressing them

import gzip
with gzip.open('strings.txt.gz', 'rt') as fp:
	print(fp.read())

# In the code above, 'rt' means to read as text (the usual thing to do)

# If you want to write named files, see 10_files_extra.py
# Instead of writing named files, use print and redirect output in the shell
# 	python3 your_program.py > your_output   # saves to a file
#	python3 your_program.py | other_program # pipes to a program


