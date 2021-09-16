#!/usr/bin/env python
import utils
import gzip
import argparse
from pysam import TabixFile
import numpy as np
import glob

def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-r',
                        dest='adj_score_dir',
                        required=True,
                        help='Path to directory containing adjusted score files')

    args = parser.parse_args()

    return args

def main():
    args = get_args()

    regions = []
    samples = []

    first_file = True

    file_i = 0

    for adj_file in glob.glob(args.adj_score_dir + '*.adj_z.bed.gz'):
        sample = adj_file.split('/')[-1].split('.')[0]
        samples.append(sample)
        with gzip.open(adj_file,'rt') as f:
            line_i = 0
            for l in f:
                A = l.rstrip().split('\t')
                region = utils.Interval(chrom=A[0],
                                        start=int(A[1]),
                                        end=int(A[2]),
                                        data=A[3:] if len(A) > 3 else None)

                if region.start == region.end:
                    continue

                interval = utils.Interval(chrom=A[0],
                                          start=int(A[1]),
                                          end=int(A[2]),
                                          data=A[3:] if len(A) > 3 else None)

                if first_file:
                    regions.append([interval])
                else:
                    regions[line_i].append(interval)

                line_i += 1
        first_file = False

        file_i += 1


    print('\t'.join(['#CHROM', 'START', 'END', 'STRAND'] + samples))

    for region in regions:
        adjs = []
        for interval in region:
            adjs.append(float(interval.data[3]))

        print('\t'.join([str(x) for x in [region[0].chrom,
                                          region[0].start,
                                          region[0].end,
                                          region[0].data[0]] +
                                          adjs]))


if __name__ == '__main__': main()
