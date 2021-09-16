#!/usr/bin/env python
import utils
import gzip
import argparse
from pysam import TabixFile
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-m',
                        dest='mosdepth_file',
                        required=True,
                        help='Path file containing target bam coverage')

    parser.add_argument('--exons_file',
                        dest='exons_file',
                        required=True,
                        help='Exon coordinate bed file')

    args = parser.parse_args()

    return args

def main():
    args = get_args()

    regions = []
    
    with gzip.open(args.exons_file,'rt') as f:
        for l in f:
            A = l.rstrip().split('\t')
            region = utils.Interval(chrom=A[0],
                                    start=int(A[1]),
                                    end=int(A[2]),
                                    data=A[3:] if len(A) > 3 else None)

            if region.start == region.end:
                continue

            regions.append(utils.Interval(chrom=A[0],
                                          start=int(A[1]),
                                          end=int(A[2]),
                                          data=A[3:] if len(A) > 3 else None))



    rates = []
    tbx = TabixFile(args.mosdepth_file)
    for region in regions:
        I = utils.get_intervals_in_region(region, args.mosdepth_file, tbx=tbx)
        total_depth = 0
        for i in I:
            total_depth += int(i.data[0])

        rate = float(total_depth)/float(region.end - region.start)
        rates.append(rate)

    #mean = np.mean(rates)
    #stdev = np.std(rates)
    
    for i in range(len(regions)):
        region = regions[i]
        rate = rates[i]
        print('\t'.join([str(x) for x in [region.chrom,
                                          region.start,
                                          region.end,
                                          region.data[0],
                                          rate]]))


if __name__ == '__main__': main()
