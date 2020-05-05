#!/usr/bin/env python3

import random

# Move the triple quotes downward to uncover each segment of code


# -------------------- FULL FUNCTION SYNTAX --------------------

# Functions often have a fixed number of positional arguments
# For example, fun1() has 2, and the order matters

def fun1(name, age):
	print(name, age)

# You call such functions with the correct order and number of parameters
fun1('Ian', 53)

# But Python also lets you call them in any order if you name them
fun1(age=53, name='Ian')

# You can also provide default arguments

def fun2(name='Ian', age=53):
	print(name, age)

fun2()
fun2(age=50)
fun2(name='Rory')
fun2('Elite', 1337)
# fun2('foo', 1, name='bar') # but this is an error

# You can also specify a function with a variable number arguments
# In which case the arguments are passed in as a tuple

def fun2(*param):  # takes 1 positional argument but one was given
	print(param)

fun2(1, 2, 3, "four")

# You can even have variable key=value arguments (dictionary)
# But just because you can, doesn't mean you should

def fun3(**param):
	print(param)

fun3(name='Ian', age=53)

# -------------------- GENERATORS --------------------

# There are times when you want to generate a list of data
# Here's a function that returns a list of random sequences

def randseq(count, length):
	seqs = []
	for i in range(count):
		chars = []
		for j in range(length):
			r = random.randint(0,3)
			if   r == 0: chars.append('A')
			elif r == 1: chars.append('C')
			elif r == 2: chars.append('G')
			else:        chars.append('T')
		seqs.append(''.join(chars))
	return seqs

# A list of a few sequences isn't a big deal

seqs = randseq(10, 5)
print(seqs)

# But imagine we want a function to create many large sequences

#seqs = randseq(1000000, 300) # size of human genome in 1M chunks

# Holding all of those sequences in memory is a lot of work
# Instead of creating all sequences up-front, generate them as needed
# Any function with a 'yield' is a generator
# A 'yield' temporarily stops execution

def seqgen(count, length):
	for i in range(count):
		chars = []
		for j in range(length):
			r = random.randint(0,3)
			if   r == 0: chars.append('A')
			elif r == 1: chars.append('C')
			elif r == 2: chars.append('G')
			else:        chars.append('T')
		yield ''.join(chars)  # create a list , iterate the list and throw it away

for seq in seqgen(10, 5):
	print(seq)

# Generators improve performance by reducing memory

# In addition to generator functions, there are generator expressions
# These look a lot like list comprehension
# Except with parentheses instead of square brackets

for r in (random.random() for i in range(10)):
	print(r)
