#!/usr/bin/env python
import pandas as pd
import sys

"""
1 - all_samples.aggregate_probe_allele_counts.txt
2 - probes file
"""

al = pd.read_csv(sys.argv[1],sep='\t',header=None)
bed = pd.read_csv(sys.argv[2],sep='\t',header=None)
samples = list( str(x) for x in al.iloc[:,0])
print('\t'.join(['chr','start','end','exon'] + samples))
for i in range(bed.shape[0]):
    print('\t'.join(list( str(x) for x in bed.iloc[i,0:4]) + [ str(x) for x in al.iloc[:,i+1]]) ) 
