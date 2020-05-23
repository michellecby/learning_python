#!/usr/bin/env python3

# Write a program that compares two files of names to find:
#	Names unique to file 1
#	Names unique to file 2
#	Names shared in both files

import sys
import biotools

file1 = sys.argv[1]
file2 = sys.argv[2]


def mkdictionary(filename):
	names = {}
	with open(filename) as fp:
		for name in fp.readlines():
			name = name.rstrip()
			names[name] = True
	return names

d1 = mkdictionary(file1)
d2 = mkdictionary(file2)

u1 = []
u2 = []
both = []

for word in d1:
	if word in d2:
		both.append(word)
	else:
		u1.append(word)

for word in d2:
	if word not in d1:
		u2.append(word)
print("Names unique to file 1:")
print(u1)
print("Names unique to file 2:")
print(u2)
print("Names shared in both files:")
print(both)

"""
python3 checklist.py file1.txt file2.txt
"""
