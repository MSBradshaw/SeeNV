import pandas as pd
import sys
import numpy as np

def reciprocal_ol(loc1, loc2, r):
    """
    @param loc1: tuple of genomic start and end cooridantes e.g. (1000,2000)
    @param loc2: tuple of genomic start and end cooridantes e.g. (1000,2000)
    @param r: reciprocal overlap threshold as a float
    return: True/False if the two locations reciprocally overlap each other
    """
    ol = min(loc1[-1], loc2[-1]) -  max(loc1[0], loc2[0])
    return  (ol / (loc1[1] - loc1[0])) >= r and (ol / (loc2[1] - loc2[0])) >= r

"""
Take in 3 bedpe files, check if there is overlap between all three breakpoints
savvy='work/{sample}.savvy.{cnv_type}.filtered.bed',
gatk='work/{sample}.gatk.{cnv_type}.filtered.bed',
cnvkit='work/{sample}.cnvkit.{cnv_type}.filtered.bed'
9       68728000        69068000        DUP     NA07347 Savvy
Inputs
1. bed file of ALL calls
2. buffer_size
3. output file
"""
try:
    bed1 = pd.read_csv(sys.argv[1], sep='\t', header=None)
except pd.errors.EmptyDataError:
    print('One of the input files is empty')
    exit(0)

buffer_size = int(sys.argv[2])
outfile = sys.argv[3]

beds = bed1
beds.columns = ['chrom','start','end','type','sample','caller']

# make sure types are all consistent
beds['chrom'] = beds['chrom'].astype(str)
beds['start'] = beds['start'].astype(int)
beds['end'] = beds['end'].astype(int)

triple_ols = []
double_ols = []
# for each chrom
for c in beds['chrom'].unique():
    sub = beds[beds['chrom'] == c]
    # compute look up matrix, m
    # n by n matrix of -1, fill each cell with the reciprocal overlap of the two calls
    m = np.full((sub.shape[0],sub.shape[0]), -1)
    for i in range(sub.shape[0]):
        for j in range(sub.shape[0]):
            # skip if it is self of has already been computed
            if i == j or m[i,j] != -1:
                continue
            # start and end points for the calls for reciprocal overlap
            i_start = int(sub.iloc[i,:]['start'])
            i_end = int(sub.iloc[i,:]['end'])
            j_start = int(sub.iloc[j,:]['start'])
            j_end = int(sub.iloc[j,:]['end'])

            # start and end of the first buffer region
            i_start1 = i_start - buffer_size
            i_end1 = i_start + buffer_size
            j_start1 = j_start - buffer_size
            j_end1 = j_start + buffer_size

            # start and end for the second buffer region
            i_start2 = i_end - buffer_size
            i_end2 = i_end + buffer_size
            j_start2 = j_end - buffer_size
            j_end2 = j_end + buffer_size

            is_reciprocal = reciprocal_ol([i_start,i_end],[j_start,j_end], r=.5)
            # these recips with r = 0.0 are for overlapping the buffers, where if there is any overlap between both buffers it is a match
            is_buffer1 = reciprocal_ol((i_start1,i_end1), (j_start1,j_end1), r=0.0)
            is_buffer2 = reciprocal_ol((i_start2,i_end2), (j_start2,j_end2), r=0.0) 
            if is_buffer1 and is_buffer2 and is_reciprocal:
                print()
                double_ols.append([])
                double_ols[-1].append('\t'.join([str(x) for x in list(sub.iloc[i,:])]))
                double_ols[-1].append('\t'.join([str(x) for x in list(sub.iloc[j,:])]))
                m[i,j] = 1
                m[j,i] = 1

lines = []
for p in double_ols:
    p.sort()
    lines.append('\t'.join(p) + '\n')
lines = set(lines)

with open(outfile,'w') as out:
    for l in lines:
        out.write(l)
