#!/usr/bin/env python3

import math

# Move the triple quotes downward to uncover each segment of code



# Floating point math is approximate math
# Because of this, you never want to use floats in a conditional statment

v = 0.1 + 0.1 + 0.1
if v == 0.3: print('equal')
else:        print('not equal', 0.3, v)

# Instead, check if numbers are close
if math.isclose(v, 0.3): print('close enough')

# It's a good idea to use assert() statements to ensure data bounds

p = 0.1 # some probability, try changing this value outside of 0-1
assert(p >= 0 and p <= 1)
