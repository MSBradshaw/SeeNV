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
                        dest='rate_dir',
                        required=True,
                        help='Path to directory containing rate files')

    args = parser.parse_args()

    return args

def main():
    args = get_args()

    regions = []

    first_file = True

    file_i = 0

    for rate_file in glob.glob(args.rate_dir + '*.probe.*.bed.gz'):
        with gzip.open(rate_file,'rt') as f:
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

    for region in regions:
        depths = []
        for interval in region:
            depths.append(float(interval.data[1]))

        print('\t'.join([str(x) for x in [region[0].chrom,
                                          region[0].start,
                                          region[0].end,
                                          region[0].data[0],
                                          np.mean(depths),
                                          np.std(depths)]]))


if __name__ == '__main__': main()
