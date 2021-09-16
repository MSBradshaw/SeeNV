#!/usr/bin/env python
import sys

"""
arguments:
1. a bedfile
2. file with a single line formated like: ['WES100', 'WES52' ...]
3. the name of an output file
"""

ids = None
for line in open(sys.argv[2],'r'):
    ids = line.replace(' ','').replace('[','').replace(']','').split(',')
    break

with open(sys.argv[3],'w') as outfile:
    for line in open(sys.argv[1],'r'):
        if sum(x in line for x in ids) > 0:
            outfile.write(line)
        
    
