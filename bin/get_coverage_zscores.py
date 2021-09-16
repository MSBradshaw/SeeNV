#!/usr/bin/env python
import utils
import gzip
import argparse
from pysam import TabixFile
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-r',
                        dest='rate_file',
                        required=True,
                        help='Path file containing coverage rate')

    parser.add_argument('-s',
                        dest='stats_file',
                        required=True,
                        help='Path file containing region stats')

    args = parser.parse_args()

    return args

def main():
    args = get_args()

    stats = []
    with gzip.open(args.stats_file,'rt') as f:
        for l in f:
            A = l.rstrip().split('\t')
            region = utils.Interval(chrom=A[0],
                                    start=int(A[1]),
                                    end=int(A[2]),
                                    data={'gene':A[3],
                                          'mean':float(A[4]),
                                          'stdev':float(A[5])})
            if region.start == region.end:
                continue

            stats.append(region)

    regions = []

    Z = []
    
    with gzip.open(args.rate_file,'rt') as f:
        i = 0
        for l in f:
            A = l.rstrip().split('\t')
            region = utils.Interval(chrom=A[0],
                                    start=int(A[1]),
                                    end=int(A[2]),
                                    data={'gene':A[3],
                                          'rate':float(A[4])})


            if region.start == region.end:
                continue

            assert region.chrom == stats[i].chrom
            assert region.start == stats[i].start
            assert region.end == stats[i].end

            z = 0
            if stats[i].data['stdev'] > 0:
                z = (region.data['rate'] - stats[i].data['mean']) / stats[i].data['stdev'] 

            region.data['z'] = z

            Z.append(z)

            regions.append(region)

            i+=1

    z_mean = np.mean(Z)
    z_stdev = np.std(Z)


    for region in regions:
        adj_z =  0
        if z_stdev > 0 :
            adj_z = (region.data['z'] - z_mean) / z_stdev
        print('\t'.join([str(x) for x in [region.chrom,
                                          region.start,
                                          region.end,
                                          region.data['gene'],
                                          region.data['rate'],
                                          region.data['z'],
                                          adj_z]]))



if __name__ == '__main__': main()
