Lesson 2
========

## Learning Objectives ##

* Variables - strings, ints, floats, None, type()
* Math - math operators and functions
* Text - string operators and functions
* Loops - iteratation

## Tutorials ##

Start with `01_variables.py`. The first thing you will learn is about
comments. That is, text you put into programs that don't actually do
anything other than provide information. Comments start with a `#` sign,
which is called number sign, pound, or hash, but not hashtag. Comments
are used during programming for debugging. One way to comment-out a
large section of code is with triple quotes.

In `01_variables.py` you will learn that variables are containers that
can hold different kinds of things, like numbers and text. Each variable
has a *type* which can be seen with the `type()` function. Variables of
different types don't always play well together, so you may need to
convert them with functions like `int()`, `float()`, or `str()`. You are
also introduced to the `None` type, which turns out to be a really
useful for debugging later.

Next, follow `02_math.py`. Here, you will learn that the typical
mathematics operations you are already familar (addition, subtraction,
multiplication, division, exponentiation) are also in Python. However,
the symbol for multiplication is `*` and exponentiation is `**`. The
standard rules for precedence are use: multiplication and division
before addition and subtraction. But what about exponentiation and other
operations like `%` (modulo)? Use parentheses when in doubt. Heck, use
parentheses even when not in doubt.

In `02_math.py` you are first introduced to the concept of libraries
(also called modules). For example, the `log2()` function is in the
`math` library. You must first `import` that math libray and then the
`log2()` function becomes available.

	import math()
	print(math.log2(1))

Next, follow `03_text.py`. Here, you will learn that some of the math
operators like `+` and `*` also work on strings. You will also see the
`len()` function for the first time, which returns the length of a
string (also the length of a list, but that comes later). You are also
introduced to the string slice syntax, which allows you to retrieve part
of a string.

Next up is `04_loops.py`. Much of programming involves looping, which is
more formally called *iteration*. A `for` loop can iterate over a range
of numbers with the `range()` function or over parts of a sequence, such
as the letters of a string. These are two fundamentally different ways
of thinking about loops, either as numbers, or temporary variables.

	dna = 'ACGT'
	
	for i in range(len(dna)):
		print(i, dna[i]) # indexing dna by number
			
	for nt in dna:
		print(nt) # temporary variable

Ultimately, you have to get really comfortable with either type of `for`
loop construction, because both are used frequently.

## Coding Exercises ##

For programming practice, try writing the following programs. The
instructions are inside each file.

1. `sumfac.py`
2. `codons.py`
3. `aa_pairs.py`
4. `frame.py`

## Python Manifest ##

### Operators

	=
	+
	-
	*
	/
	**
	%
	//
	\
	
### Keywords

	for
	in
	
### Functions

	float()
	int()
	len()
	range()
	str()
	type()
	
### Libraries

	math.e
	math.pi
	math.inf
	math.nan
	math.tau
	math.ceil()
	math.factorial()
	math.log()
	math.log2()
	math.floor()
	math.sqrt()

## File Manifest ##

	00README.md
	01_variables.py
	02_math.py
	03_text.py
	04_loops.py
	aa_pairs.py
	codons.py
	frame.py
	sumfac.py
