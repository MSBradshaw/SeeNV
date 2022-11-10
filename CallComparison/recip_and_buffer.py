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
1. bed file
2. bed file
3. bed file
4. buffer_size
5. output file
6. double overlap output file
"""
try:
    bed1 = pd.read_csv(sys.argv[1], sep='\t', header=None)
    bed2 = pd.read_csv(sys.argv[2], sep='\t', header=None)
    bed3 = pd.read_csv(sys.argv[3], sep='\t', header=None)
except pd.errors.EmptyDataError:
    print('One of the input files is empty')
    exit(0)

buffer_size = int(sys.argv[4])
outfile = sys.argv[5]
outfile_doubles = sys.argv[6]

beds = pd.concat([bed1,bed2,bed3])
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
                m[i,j] = 1
                m[j,i] = 1
            else:
                m[i,j] = 0
                m[j,i] = 0
    # find all sets of 3 calls
    for i in range(sub.shape[0]):
        for j in range(i,sub.shape[0]):
            # if the first two don't match, skip
            if m[i,j] != 1:
                continue
            double_ols.append([])
            double_ols[-1].append('\t'.join([str(x) for x in list(sub.iloc[i,:])]))
            double_ols[-1].append('\t'.join([str(x) for x in list(sub.iloc[j,:])]))
            for k in range(j,sub.shape[0]):
                # if the last two match then they all match
                if m[j,k] == 1:
                    triple_ols.append([])
                    triple_ols[-1].append('\t'.join([str(x) for x in list(sub.iloc[i,:])]))
                    triple_ols[-1].append('\t'.join([str(x) for x in list(sub.iloc[j,:])]))
                    triple_ols[-1].append('\t'.join([str(x) for x in list(sub.iloc[k,:])]))
lines = []
for pairs in triple_ols:
    pairs.sort()
    lines.append('\t'.join(pairs) + '\n')
lines = set(lines)

double_lines = []
for p in double_ols:
    p.sort()
    double_lines.append('\t'.join(p) + '\n')
double_lines = set(double_lines)

with open(outfile,'w') as out:
    for l in lines:
        out.write(l)

with open(outfile_doubles,'w') as out:
    for l in double_lines:
        out.write(l)
