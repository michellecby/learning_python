#!/usr/bin/env python3

# Move the triple quotes downward to uncover each segment of code



# -------------------- ZIP FUNCTION --------------------

# The zip() function lets you loop through iterables in parallel
# The zip() ends when the first container ends

names = ('Nigel', 'David', 'Derek') # tuple
ages = [52, 53, 49, 1, 2, 3]        # list

for (name, age) in zip(names, ages):
	print(name, age)


# -------------------- LIST COMPREHENSION --------------------

# Consider the following initialization code:

data = []
for i in range(10):
	data.append(0)
print(data)

# You can write this more succinctly with the * operator

data = [0] * 10
print(data)

# Another alternative is list comprehension

data = [0 for i in range(10)]
print(data)

# List comprehension gets even more interesting...
# Here's a slightly more complex initialization

squares = []
for i in range(10):
	squares.append(i ** 2)
print(squares)

# List comprehension turns 3 lines into 1

squares = [i ** 2 for i in range(10)]
print(squares)

# You can also include a conditional
# First the longer syntax

square3 = []
for i in range(10):
	if i % 3 == 0:
		square3.append(i ** 2)
print(square3)

# Now list comprehension

square3 = [i ** 2 for i in range(10) if i % 3 == 0]
print(square3)
