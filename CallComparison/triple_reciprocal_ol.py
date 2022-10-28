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
savvy='work/{sample}.savvy.{cnv_type}.filtered.bed',
gatk='work/{sample}.gatk.{cnv_type}.filtered.bed',
cnvkit='work/{sample}.cnvkit.{cnv_type}.filtered.bed'
9       68728000        69068000        DUP     NA07347 Savvy
Inputs
1. bed file
2. bed file
3. bed file
4. output file
"""

bed1 = pd.read_csv(sys.argv[1], sep='\t', header=None)
bed2 = pd.read_csv(sys.argv[2], sep='\t', header=None)
bed3 = pd.read_csv(sys.argv[3], sep='\t', header=None)
outfile = sys.argv[4]

beds = pd.concat([bed1,bed2,bed3])
beds.columns = ['chrom','start','end','type','sample','caller']
triple_ols = []
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
            # calc overlap between the two calls
            is_reciprocal = reciprocal_ol((int(sub.iloc[i,:]['start']),int(sub.iloc[i,:]['end'])), (int(sub.iloc[j,:]['start']),int(sub.iloc[j,:]['end'])), r=.5)
            if is_reciprocal:
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
