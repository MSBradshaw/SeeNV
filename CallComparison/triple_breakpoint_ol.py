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

beds = pd.concat([bed1,bed2,bed3])
beds.columns = ['chrom1','start1','end1','chrom2','start2','end2','sample','score','strand1','strand2','call_type','caller']
triple_ols = []
# for each chrom
for c in beds['chrom1'].unique():
    sub = beds[beds['chrom1'] == c]
    # compute look up matrix, m
    # n by n matrix of -1, fill each cell with the reciprocal overlap of the two calls
    m = np.full((sub.shape[0],sub.shape[0]), -1)
    for i in range(sub.shape[0]):
        for j in range(sub.shape[0]):
            # skip if it is self of has already been computed
            if i == j or m[i,j] != -1:
                continue
            # calc overlap between the two calls with buffer overlap, r=0.0 means any overlap will count
            i_start1 = int(sub.iloc[i,:]['start1']) - buffer_size
            i_end1 = int(sub.iloc[i,:]['end1']) + buffer_size
            j_start1 = int(sub.iloc[j,:]['start1']) - buffer_size
            j_end1 = int(sub.iloc[j,:]['end1']) + buffer_size

            i_start2 = int(sub.iloc[i,:]['start2']) - buffer_size
            i_end2 = int(sub.iloc[i,:]['end2']) + buffer_size
            j_start2 = int(sub.iloc[j,:]['start2']) - buffer_size
            j_end2 = int(sub.iloc[j,:]['end2']) + buffer_size

            is_reciprocal1 = reciprocal_ol((i_start1,i_end1), (j_start1,j_end1), r=0.0)
            is_reciprocal2 = reciprocal_ol((i_start2,i_end2), (j_start2,j_end2), r=0.0) 
            if is_reciprocal1 and is_reciprocal2:
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

with open(outfile,'w') as out:
    for l in lines:
        out.write(l)
