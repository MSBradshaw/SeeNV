import sys

"""
1   16877000    16894000    DUP HG00124 Savvy   1   16888814    16940057    NBPF1   .   -   5186
10  37441000    37482000    DUP HG00124 Savvy   10  37414785    37673039    ANKRD30A    .   +   41000
11  55370000    55420000    DUP HG00124 Savvy   11  55370830    55371874    OR4C11  .   -   1044
"""
res = {}
for line in sys.stdin:
    row = line.strip().split('\t')
    uid = '\t'.join(row[0:6])
    if uid not in res:
        res[uid] = []
    # call location -> [gene1,gene2,gene3]
    res[uid].append(row[9])

for uid in res.keys():
    print(uid + '\t' + ','.join(res[uid]))
