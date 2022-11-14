import sys

# 1 12853000    12944000    DEL HG00124 Savvy   PRAMEF1,PRAMEF11,HNRNPCL1,PRAMEF2,PRAMEF4

# gene: [uid,uid]
res = {}
for line in sys.stdin:
    row = line.strip().split('\t')
    uid = '\t'.join(row[0:6])
    genes = row[6].split(',')
    for g in genes:
        if g not in res:
            res[g] = []
            res[g].append(uid)
        else:
            res[g].append(uid)
            

# count callers in uid at each gene
# how many calls affect each gene
# intentionally leaving a blank so 3 layers of loops will do all possible combos
callers = ['GATK','CNVkit','Savvy','\t', '\t']
counts = {}
for g in res.keys():
    #print(g,res[g])
    for i in range(len(callers)):
        for j in range(i,len(callers)):
            if i == j:
                continue
            for k in range(j,len(callers)):
                if j == k or i == k:
                    continue
                ci = callers[i]
                cj = callers[j]
                ck = callers[k]
                if sum([ci in x for x in res[g]]) > 0 and sum([cj in x for x in res[g]]) > 0 and sum([ck in x for x in res[g]]) > 0:
                    caller_combo = [ci,cj,ck]
                    caller_combo.sort()
                    if str(caller_combo) not in counts:
                        counts[str(caller_combo)] = 0
                    counts[str(caller_combo)] += 1

for combo in counts.keys():
    print(combo, counts[combo])

