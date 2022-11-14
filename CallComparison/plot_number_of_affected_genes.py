import sys
import matplotlib.pyplot as plt

"""
1. output file name
2. output file 2 for zoomed in x axis
"""

# cat work/GeneOverlap/*.DEL.filtered.grouped.bed
# 8 11970074    12291651    DUP HG00124 CNVkit  ZNF705D,USP17L2,FAM85A,FAM66D,FAM86B1,DEFB130,FAM86B2

counts = []
for line in sys.stdin:
    row = line.strip().split('\t')
    counts.append(len(row[6].split(',')))
    if counts[-1] > 1000:
        print(row[:6],counts[-1])

plt.hist(counts,bins=20)
plt.yscale('log')
plt.ylabel('count')
plt.xlabel('# genes affected')
plt.savefig(sys.argv[1])
plt.clf()

small_counts = [x for x in counts if x < 200]
plt.hist(small_counts,bins=100)
plt.yscale('log')
plt.ylabel('count')
plt.xlabel('# genes affected')
plt.savefig(sys.argv[2])
