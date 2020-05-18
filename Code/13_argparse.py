#!/usr/bin/env python3

import argparse

# This is a template you can use for processing command line arguments
# Change the names of the arguments to something descriptive

# Commandline: python3 13_argparse.py --rstr yo --rint 1 --rfloat 0.1
# Also try:
#	python3 13_argparse.py
#	python3 13_argparse.py --help



# setup
parser = argparse.ArgumentParser(
	description='Brief description of program.')
# required arguments
parser.add_argument('--rstr', required=True, type=str,
	metavar='<str>', help='required string argument')
parser.add_argument('--rint', required=True, type=int,
	metavar='<int>', help='required integer argument')
parser.add_argument('--rfloat', required=True, type=float,
	metavar='<float>', help='required floating point argument')
# optional arguments with default parameters
parser.add_argument('--dstr', required=False, type=str, default='hello',
	metavar='<str>', help='optional string argument [%(default)s]')
parser.add_argument('--dint', required=False, type=int, default=1,
	metavar='<int>', help='optional integer argument [%(default)i]')
parser.add_argument('--dfloat', required=False, type=float, default=3.14,
	metavar='<float>', help='optional floating point argument [%(default)f]')
# switches
parser.add_argument('--switch', action='store_true',
	help='on/off switch')
# finalization
arg = parser.parse_args()

# testing
print(arg.rstr, arg.rint, arg.rfloat)
print(arg.dstr, arg.dint, arg.dfloat)
if arg.switch: print('switch on')
else:          print('switch off')
