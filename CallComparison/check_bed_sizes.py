import sys

small_ones = 0
for line in open(sys.argv[1],'r'):
    row = line.strip().split('\t')
    if int(row[2]) - int(row[1]) < 1000:
        small_ones += 1
print(small_ones)
