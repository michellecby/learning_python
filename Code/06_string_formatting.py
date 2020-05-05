#!/usr/bin/env python3

# Move the triple quotes downward to uncover each segment of code

"""

txt = 'Ian'
num = 3/11

# Previously, we have used the print() function with commas

print(txt, num)

# What if we want to control the way the text looks?
# For example, what if we want exactly 3 decimal places?
# There are 3 distinct ways to format strings in python

# Method 1: printf()
# printf() is the name of an old C function
# The syntax is well-known among programmers but looks odd to novices

print('%s %.3f' % (txt, num))                # %s string, %f float
print('%s %.3f %d %e' % (txt, num, 2.1, .1)) # %d integer, %e scientific

# Method 2: str.format()
# Strings are objects with built-in methods (i.e. functions)
# upper() and lower() are some simple examples of string methods

print(txt.upper(), txt.lower())

# The format() method is a powerful way to control string formatting

print('{} {}'.format(txt, num))
print('{} {:.3f}'.format(txt, num))

# Method 3: f-strings
# f-strings let you interpolate variables inside curly brackets

print(f'{txt} {num}')
print(f'{txt} {num:.3f}')

# You can even interpolate python code

print(f'{2+2} {1/7:.5f} {len(txt)}')

# The examples here are but the tip of a very large iceberg
# Each formatting method has many more options
# Check documentation online for more information

"""
