#!/usr/bin/env python
import collections
from pysam import VariantFile
from pysam import AlignmentFile
from pysam import TabixFile
import gzip
import numpy as np

Interval = collections.namedtuple('Interval', 'chrom start end data')

def get_gene_cord_map(gene_bed):
    genes = {}

    f =  gzip.open(gene_bed, 'rb')
    for l in f:
        A = l.rstrip().split()
        cord = Interval(chrom=A[0].decode('UTF-8'),
                        start=int(A[1]),
                        end=int(A[2]),
                        data=None)
        name = A[3].decode('UTF-8')
        genes[name]=cord

    f.close()

    return genes

def get_header(bed_file, tbx=None):
    if tbx is None:
        tbx = TabixFile(bed_file)
    return tbx.header

def get_intervals_in_region(target_interval, bed_file, tbx=None, ignore=None):
    if tbx is None:
        tbx = TabixFile(bed_file)
    intervals = []
    for row in tbx.fetch(target_interval.chrom,
                         target_interval.start,
                         target_interval.end):
        if ignore is not None and ignore in row: continue
        A = row.rstrip().split()
        cord = Interval(chrom=A[0],
                        start=int(A[1]),
                        end=int(A[2]),
                        data=A[3:] if len(A)>3 else None)
        intervals.append(cord)
    return intervals

def get_gene_exons(target_gene, gene_bed_file, exons_bed_file):
    gene_cord_map = get_gene_cord_map(gene_bed_file)
    gene_cord = gene_cord_map[target_gene]
    exons = get_intervals_in_region(gene_cord, exons_bed_file)

    return exons

def get_gt_counts(sample, region, vcf_file=None, vcf_handle=None):

    vcf = None
    
    if vcf_handle is None:
        vcf = VariantFile(vcf_file)
    else:
        vcf = vcf_handle

    gt_count = {'HOM_REF':0,
                'HET':0,
                'HOM_ALT':0,
                'UNK':0}

    for rec in vcf.fetch(region.chrom, region.start, region.end):
        gt = rec.samples[sample]['GT']

        if gt[0] is None or gt[1] is None:
            gt_count['UNK'] += 1
        elif gt[0] == gt[1]:
            if gt[0] == 0:
                gt_count['HOM_REF'] += 1
            else:
                gt_count['HOM_ALT'] += 1
        else:
            gt_count['HET'] += 1

    return gt_count

def get_gt_stats(region, target_sample=None, vcf_file=None, vcf_handle=None):

    vcf = None
    
    if vcf_handle is None:
        vcf = VariantFile(vcf_file)
    else:
        vcf = vcf_handle

    samples = None
    target_counts = [0,0,0]

    for rec in vcf.fetch(region.chrom, region.start, region.end):

        if samples is None:
            samples = {}
            for sample in rec.samples:
                samples[sample] = [0,0,0]

        for sample in rec.samples:
            gt = rec.samples[sample]['GT']
            if gt[0] == None:
                samples[sample][2] += 1
            elif gt[0] == gt[1]:
                samples[sample][0] += 1
            else:
                samples[sample][1] += 1

        gt = rec.samples[target_sample]['GT']
        if gt[0] == None:
            target_counts[2] += 1
        elif gt[0] == gt[1]:
            target_counts[0] += 1
        else:
            target_counts[1] += 1

    if samples is None:
        return None

    counts = [[],[],[]]

    for sample in samples:
        for i in range(3):
            counts[i].append(samples[sample][i])

    count_stats = []

    for i in range(3):
        m = np.mean(counts[i])
        s = np.std(counts[i])
        x = target_counts[i]
        z = 0
        if s > 0 : z = (x-m)/s
        count_stats.append( [z,x,m,s] )


    return count_stats



#    r = []
#    gts = []
#    ads = []
#    dps = []
#    gqs = [] 
#
#
#    S = None
#
#
#    for rec in vcf.fetch(region.chrom, region.start, region.end):
#        samples = []
#        if target_sample:
#            samples = [target_sample]
#        else:
#            samples = rec.sample
#
#        for sample in samples:
#            ad = rec.samples[sample]['AD']
#            gt = rec.samples[sample]['GT']
#            dp = rec.samples[sample]['DP']
#            gq = rec.samples[sample]['GQ']
#            if gt[0] and gt[1]:
#                gts.append(gt)
#                ads.append(ad)
#                dps.append(dp)
#                gqs.append(gq)
#    return (gts, ads, dps, gqs)

def get_coverage(region, bam_file=None, bam_handle=None):

    bam = None
    
    if bam_handle is None:
        bam = AlignmentFile(bam_file, 'rb')
    else:
        bam = bam_handle

    iter = bam.fetch(region.chrom, region.start, region.end)
    count = 0
    for x in iter:
        count += 1

    return count
