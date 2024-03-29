#!/usr/bin/env python
import utils
import gzip
import argparse
from pysam import TabixFile
import numpy as np
import warnings

def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-m',
                        dest='mosdepth_file',
                        required=True,
                        help='Path file containing target bam coverage')

    parser.add_argument('--num_reads',
                        dest='num_reads',
                        required=True,
                        help='Number of total reads in BAM/CRAM file')

    parser.add_argument('--regions_file',
                        dest='regions_file',
                        required=True,
                        help='Region coordinate bed file')

    args = parser.parse_args()

    return args

def main():
    args = get_args()
    num_reads = None
    for line in open(args.num_reads,'r'):
        num_reads = float(line.strip())
    rpm_factor = num_reads/1000000

    regions = []
    
    with gzip.open(args.regions_file,'rt') as f:
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
        try:
            I = utils.get_intervals_in_region(region, args.mosdepth_file, tbx=tbx)
        except ValueError as e:
            warnings.warn('Warning, adding np.nan as rate for this region:'+str(e))
            rates.append(np.nan)
            continue
        total_depth = 0
        for i in I:
            total_depth += int(i.data[0])

        rate = float(total_depth)/rpm_factor/float(region.end-region.start)
        rates.append(rate)

    #mean = np.mean(rates)
    #stdev = np.std(rates)
    
    for i in range(len(regions)):
        region = regions[i]
        rate = rates[i]
        if region.data != None:
            gene_name = region.data[0]
        else:
            gene_name = 'No_gene_name'
        print('\t'.join([str(x) for x in [region.chrom,
                                          region.start,
                                          region.end,
                                          gene_name,
                                          rate]]))


if __name__ == '__main__': main()
