Lesson 5
========

## Learning Objectives ##

* Open and read from named files, both normal and compressed
* Search for items inside lists
* Define functions
* Iterate through fasta files

## Tutorials ##

Start with `10_files.py`, where you will learn how to open and read from
named files, for example, those named on the command line. You can also
open and read from compressed files without first decompressing them
using the `gzip` library. If you want to learn how to write or append to
files, check out `10_files_extra.py`.

The next tutorial is really short, `11_in.py`, but it has two very
important concepts: the use of `in` for searching, and the the `time`
library for benchmarking.

The main part of the lesson is `12_functions.py`. Functions are units of
code re-use and are one of the most important concepts in programming.
In addition to the instruction, you will also find a `read_fasta()`
function inside `12_functions.py`, which we will copy-paste into other
programs (for now, later we will put it in a proper library).

Python function syntax is highly customizable, but sort of confusing.
For that reason, you can find a more complete explanation in
`12_functions_extra.py`. Here also you will find *generator* functions,
which let you turn functions into iterators.

## Coding Exercises ##

For programming practice, try writing the following programs. The
instructions are inside each file.

1. `exonic_snps.py`
2. `pepsearch.py`
3. `fasta_stats.py`
4. `rand_fasta.py`
5. `longest_orf.py`
6. `transmembrane.py`
7. `longest_cds.py` (requires `12_functions_extra.py`)

## Python Manifest ##

### Keywords

	def
	return
	with
	yield

### Functions

	open()

### Methods

	file.read()
	file.readline()
	file.readlines()
	
	time.perf_counter()

### Libraries

	gzip.open()

## File Manifest ##

	00README.md
	10_files.py
	10_files_extra.py
	11_in.py
	12_functions.py
	12_functions_extra.py
	exonic_snps.py
	exons.txt.gz
	fasta_stats.py
	longest_cds.py
	longest_orf.py
	pepsearch.py
	proteins.fasta.gz
	rand_fasta.py
	snps.txt.gz
	strings.txt.gz
	transcripts.fasta.gz
	transmembrane.py

