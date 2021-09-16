#!/usr/bin/env python
import sys

"""
input: tab seporated bed format with one additional field, the gene symbol
output: tab seporated bed format followed by gene symbol and exon number
"""

current = None
count = 0
for line in sys.stdin:
    row = line.strip().split('\t')
    if row[-1] == current:
        count += 1
    else:
        current  = row[-1]
        count = 1
    print('\t'.join(row + [str(count)]))
