#!/opt/miniconda/envs/chco/bin/python
import argparse
import utils
import pandas as pd
import numpy as np
import pysam
import statistics
from random import randrange

def get_args():
    parser = argparse.ArgumentParser()


    parser.add_argument('--all_calls',
                        dest='all_calls',
                        help='file with path to BGZIPED and TABIXED Savvy calls set on each line')

    parser.add_argument('--exons',
                        dest='exons',
                        help='BGZIPED and TABIXED bed file with location of all exons')

    parser.add_argument('--sample',
                        dest='sample',
                        help='Sample name')

    parser.add_argument('--region',
                        dest='region',
                        help='Target region')


    parser.add_argument('--window',
                        dest='window',
                        type=int,
                        default=50000,
                        help='Window (default 50000)')

    parser.add_argument("--depth",
                        dest="depth",
                        help="file with depth/rpm data")

    args = parser.parse_args()
    
    return args


def get_all_savvy_calls_in_target(savvy_bed_file, target_interval, ignore=None):
    calls = utils.get_intervals_in_region(target_interval, savvy_bed_file,ignore=ignore)
    sample_calls = []
    for call in calls:
        curr_sample = call.data[6].split('.')[0]
        sample_calls.append(call)
    return(sample_calls)

def main():

    args = get_args()
    call_colors = ['#37dc94', '#8931EF', '#F2CA19', '#FF00BD', '#0057E9', '#87E911', '#E11845']

    chrom=args.region.split(':')[0]

    target_region = utils.Interval(
            chrom=chrom,
            start=max(0,
                      int(args.region.split(':')[1].split('-')[0]) - args.window),
            end=int(args.region.split(':')[1].split('-')[1]) + args.window,
            data=None)


    exons = utils.get_intervals_in_region(target_region,
                                           args.exons)

    if args.exons:
        probes = []
        for line in open(args.depth,'r'):
            row = line.strip().split('\t')
            if str(row[0]) != str(target_region.chrom):
                continue
            probes.append(utils.Interval(chrom=row[0],start=int(row[1]),end=int(row[2]),data=None))

        legend_elements = []
        colors = call_colors
        all_max_ys = []
        for i,line in enumerate(open(args.all_calls,'r')):
            f = line.strip()
            try:
                density_of_calls_at_each_probe = [ len(get_all_savvy_calls_in_target(f,x,ignore=args.sample)) for x in probes]
            except ValueError:
                print('ValueError in file ' + f)
                density_of_calls_at_each_probe = [0 for x in probes]
            x_vals = [ x.start for x in probes]
            indexes = [j for j,x in enumerate(density_of_calls_at_each_probe) if x > 0]
            x_vals = [ x_vals[j] for j in indexes]
            y_vals = [ density_of_calls_at_each_probe[i] for i in indexes]
            if len(y_vals) > 0:
                all_max_ys.append(max(y_vals))
        if len(all_max_ys) > 0:
            print(max(all_max_ys))
        else:
            print(0)

main()
