#!/usr/bin/env python
import os 
import pysam
import sys
import re



f = sys.argv[2]
if '.tbi' in f:
    f = '.'.join(f.split('.')[:-1])

if len(sys.argv) >= 4:
    sample = sys.argv[3]
else:
    sample = f.split('/')[-1].split('.')[0]


res = []
tabixfile = pysam.TabixFile(f)
for i,line in enumerate(open(sys.argv[1],'r')):
    row = line.strip().split('\t')
    try:
        a = list(tabixfile.fetch(row[0],int(row[1]),int(row[2])))
    except ValueError:
        print('Error')
        continue
    if len(a) > 0:
        alt_counts = sum( int(x.strip().split('\t')[4]) for x in a)
        alt_counts_norm = alt_counts / float(len(a))
        probe = row[3]
        res.append(alt_counts_norm)
    else:
        res.append(-1)
print('\t'.join([sample] + [str(x) for x in res]))
