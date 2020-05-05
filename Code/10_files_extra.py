#!/usr/bin/env python3

import sys

# Move the triple quotes downward to uncover each segment of code


# The 'open' function defaults to reading
"""
fp = open(sys.argv[1]) # open for reading
print(fp.readline())
fp.close()             # close file pointers after you're done

# But you can also open files under other modes
# Mode:
#	r:  reading
#	w:  writing (overwrites anything already there)
#	x:  create new file, then write (will not overwrite old file)
#	a:  append to file
#	r+: read and write same file

fp = open(sys.argv[1], 'r') # open for reading, explicit syntax
print(fp.readline())
fp.close()
"""
fp = open('for_example/new_file.txt', 'x') # open for writing, will not overwrite
fp.write('new file!\n')
fp.close()

fp = open('for_example/file1.txt', 'w') # open for writing
fp.write('whatever was here has been destroyed\n')
fp.close()

fp = open('for_example/file1.txt', 'a') # open for appending
fp.write('stuff added\n')
fp.close()

fp = open('for_example/file1.txt', 'r+') # open for reading and writing
print(fp.readline())
fp.write('more stuff added\n')
print(fp.readline())
fp.close()

# If you don't close() your file pointers, bad stuff can happen
# So use the 'with' syntax whenever possible because it closes automatically

with open('/tmp/foo', 'a') as fp:
	fp.write('with ftw\n')

# Just as you can read from gz files, you can write to them as well
# Usually you don't need to do that
# Read the official documentation if that interests you
