#!/opt/miniconda/envs/chco/bin/python
import argparse
import utils
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pysam import VariantFile
import pandas as pd
import numpy as np
import pysam
import statistics
from random import randrange
from matplotlib.lines import Line2D

#{{{ def get_args():
def get_args():
    parser = argparse.ArgumentParser()


    parser.add_argument('--savvy_calls',
                        dest='savvy_calls',
                        help='BGZIPED and TABIXED Savvy calls set')

    parser.add_argument('--all_calls',
                        dest='all_calls',
                        help='file with path to BGZIPED and TABIXED Savvy calls set on each line')

    parser.add_argument('--scores',
                        dest='scores',
                        help='BGZIPED and TABIXED merged and adjusted scores')

    parser.add_argument('--exons',
                        dest='exons',
                        help='BGZIPED and TABIXED bed file with location of all exons')

    parser.add_argument('--vcf',
                        dest='vcf',
                        help='SNP/INDEL VCF"')

    parser.add_argument('--alt_allele_counts',
                        dest='alt_allele_counts',
                        help='tsv file with a normalized alterate allele count for each sample at each probe')

    parser.add_argument('--alt_allele_headers',
                        dest='alt_allele_headers',
                        help='headers for the tsv file with a normalized alterate allele count for each sample at each probe')

#    parser.add_argument('--cnvkit_calls',
#                        dest='cnvkit_calls',
#                        help='BGZIPED and TABIXED CNVKit calls set')
#
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

    parser.add_argument('-o',
                        dest='outfile',
                        required=True,
                        help='Output file name')

    parser.add_argument('--legend_loc',
                        dest='legend_loc',
                        default='best',
                        help='Legend loctaion')
 
    parser.add_argument('--width',
                        dest='width',
                        type=float,
                        default=5,
                        help='Plot width (default 5)')

    parser.add_argument('--height',
                        dest='height',
                        type=float,
                        default=5,
                        help='Plot height (default 5)')

    parser.add_argument('--tick_line_length',
                        dest='tick_line_length',
                        type=float,
                        default=2,
                        help='Tick line width')

    parser.add_argument('--tick_line_width',
                       dest='tick_line_width',
                       type=float,
                       default=0.5,
                       help='Tick line width')

    parser.add_argument('--axis_line_width',
                       dest='axis_line_width',
                       type=float,
                       default=0.5,
                       help='Axis line width')

    parser.add_argument('--axis_label_size',
                        dest='axis_label_size',
                        type=int,
                        default=8,
                        help='Axis label font size')

    parser.add_argument('--tick_label_size',
                        dest='tick_label_size',
                        type=int,
                        default=8,
                        help='Axis tick label font size')

    parser.add_argument('--x_label',
                        dest='x_label',
                        help='X axis label')

    parser.add_argument('--y_label',
                        dest='y_label',
                        help='Y axis label')

    parser.add_argument("--title",
                        dest="title",
                        help="Plot title (title or title;size;location)")

    parser.add_argument("--depth",
                        dest="depth",
                        help="file with depth/rpm data")

    parser.add_argument("--gnomad_sv",
                        dest="gnomad_sv",
                        help="bgzipped and tabixed bed file with the gnomad sv sites")

    parser.add_argument("--additional_beds",
                        dest="additional_beds",
                        help="a basic text file containing a list of bgzipped and tabixed bed files to be included")

    parser.add_argument('--label_exons',dest='label_exons', action='store_true',help="set to label each exon with it's gene name and exon number. requires that the file in --exons has gene and exon number listed in the two fields immediatelty after teh normal bed format items")

    parser.add_argument("--max_num_calls",
                        dest="max_num_calls",
                        help="what to max the max y lim on the Num. Calls plot. Can all be a file path to a file with a single line formatted like ['1','value thing','2' ... ]")

    args = parser.parse_args()
    
    return args
#}}}

def process_max_num_calls(file_name):
    try:
        info = None
        for line in open(file_name,'r'):
            info = line
            break
        info = info.replace('[','').replace(']','').replace(' ','')
        data = info.split(',')
        max_n_calls = max([ float(x) for x in data if x.isnumeric()])
        return max_n_calls
    except FileNotFoundError:
        if file_name.isnumeric():
            return float(file_name)
        else:
            raise Exception("the parameter --max_num_calls is not a properly formated file or numerical value")
            return -1


def plot_coverage_stats(_ax, _args, _call_colors, _non_samp_pos, _non_samp_xs, _regions_pos, _means, _stdevs, _xs, _chrom):
    if _args.title:
        if ';' in _args.title:
            _text, _size, _loc = _args.title.split(';')
            _fontdict = {'fontsize':int(_size)}
            _ax.set_title(_text, fontdict=_fontdict, loc='right')
        else:
            _ax.set_title(_args.title, loc='right')
        
    _lns2b = _ax.plot(_non_samp_pos,
        _non_samp_xs,
        '-o',
        lw = 0,
        markersize=1,
        label='Other Sample Coverage',
        c='#e2e2e2', alpha=.3)
        
    _lns1 = _ax.plot(_regions_pos,
        _means,
        marker='_',
        markersize=3,
        lw = 0,
        label='Mean Pop. Coverage',
        c='#5683d6')

    for i in range(len(_means)):
        _ax.plot([_regions_pos[i],_regions_pos[i]],
            [_means[i]-_stdevs[i],_means[i]+_stdevs[i]],
            lw=1,
            color='#5683d6',
            alpha=0.5)


    _lns2 = _ax.plot(_regions_pos,
        _xs,
        '-o',
        lw = 0,
        markersize=1,
        label='Sample Coverage',
        c='#5683d6')
        
    _ax.spines['left'].set_visible(False)

    format_axis(_ax, _args)
    #_ax.set_ylabel('Z-score', fontsize=6)
    _ax.set_title('Z-score', fontsize=6, x=0.0, fontdict={'horizontalalignment': 'left'})
    _xmin,_xmax = [int(x) for x in _ax.get_xlim()]
    _target_interval = utils.Interval(chrom=_chrom, start=max(0,_xmin), end=_xmax, data=None)
    _savvy_calls = []
    if _args.savvy_calls and _args.sample:
        _xmin,_xmax = [int(x) for x in ax.get_xlim()]
        _target_interval = utils.Interval(chrom=_chrom, start=max(0,_xmin), end=_xmax, data=None)
        _savvy_calls = get_savvy_calls(_args.savvy_calls, _args.sample, _target_interval)
        _all_savvy_calls = get_all_savvy_calls_in_target(_args.savvy_calls,_target_interval)

    mark_intervals(_ax, _savvy_calls, _call_colors[0],a=0.2)

    if _args.all_calls:
        _all_calls = []
        _number_of_calls_to_plot = 0
        for _line in open(_args.all_calls,'r'):
            _line = _line.strip()
            try:
                _tbx_calls = pysam.TabixFile(_line)
                _raw_calls = _tbx_calls.fetch(_target_interval.chrom, _target_interval.start, _target_interval.end)
                _calls = []
                for l in _raw_calls:
                    if _args.sample not in l: continue
                    _row = l.strip().split('\t')
                    _calls.append(utils.Interval(chrom=_row[0],
                            start=int(_row[1]),
                            end=int(_row[2]),
                            data=None))

            except ValueError:
                print('ValueError with file ' + _line)
                _calls = []
            #if there is 1 or more calls, add 1 to the count
            if len(_calls) > 0:
                _number_of_calls_to_plot += len(_calls)
            _all_calls.append(_calls)

            _box_num = 0
            for i,_calls in enumerate(_all_calls):
                if len(_calls) == 0: continue
                for _call in _calls:
                    if _number_of_calls_to_plot == 0: break 
                    _ymin = 1.0 / _number_of_calls_to_plot * _box_num
                    _ymax = _ymin + (1.0 / _number_of_calls_to_plot)
                    _color_i = i % len(_call_colors)
                    if i > len(_call_colors):
                        print('Warning, too many samples, colors are being recycled')
                    mark_intervals(_ax, [_call], _call_colors[_color_i],a=0.2,ymin=_ymin,ymax=_ymax)
                    _box_num += 1
        _ax.set_ylim([-7, 7])
        _ax.tick_params(axis='both', which='major', labelsize=4,length=2,width=.5)
        _ax.tick_params(axis='both', which='minor', labelsize=4,length=2,width=.5)
        return _ax, _lns1, _lns2b, _lns2, _all_calls, _savvy_calls, _xmin, _xmax, _target_interval 



def get_file_name(f):
    """
    removes file endings and file path information
    example: '/my/file/path.txt' -> 'path'
    f: str, file to have
    """
    return f.split('/')[-1].split('.')[0]

def get_region_alt_counts(chr,st,fin,bed_tbx,sample_col):
    res = bed_tbx.fetch(chr,st,fin)
    res_dict = {'chr':[],'start':[],'end':[],'avg':[],'sample':[],'zscore':[],'std':[]}
    for r in res:
        row = r.strip().split('\t')
        vals = [float(x) for x in row[4:]]
        mean = sum( vals) / len(vals)
        std = statistics.stdev(vals)
        if std != 0:
            zscore = (float(row[sample_col]) - mean) / std
        else:
            zscore = None
        res_dict['chr'].append(row[0])
        res_dict['start'].append(row[1])
        res_dict['end'].append(row[2])
        res_dict['avg'].append(mean)
        res_dict['sample'].append(row[sample_col])
        res_dict['zscore'].append(zscore)
        res_dict['std'].append(std)
    return pd.DataFrame(res_dict)

def format_axis(ax, args):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['bottom'].set_linewidth(args.axis_line_width)
    ax.spines['left'].set_linewidth(args.axis_line_width)
    ax.spines['right'].set_linewidth(args.axis_line_width)
    ax.tick_params(axis='both',
                   which='major',
                   labelsize=args.axis_label_size,
                   width=args.tick_line_width,
                   length=args.tick_line_length)

def mark_intervals(ax, intervals, color,a=0.25,ymin=0,ymax=1):
    for interval in intervals:
        ax.axvspan(max(ax.get_xlim()[0], interval.start),
                   min(ax.get_xlim()[1], interval.end),
                   ymin,ymax,
                   alpha=a, color=color)

def label_exons(ax, intervals):
    texts = []
    for interval in intervals:
        texts.append(ax.text(interval.start, 1, interval.data[0] + ' ' + interval.data[1], fontsize=3,rotation=90))
    return texts

def get_savvy_calls(savvy_bed_file, target_sample, target_interval):
    calls = utils.get_intervals_in_region(target_interval, savvy_bed_file)
    sample_calls = []
    for call in calls:
        curr_sample = call.data[6].split('.')[0]
        if curr_sample == target_sample:
            sample_calls.append(call)
    return(sample_calls)

def get_all_savvy_calls_in_target(savvy_bed_file, target_interval, ignore=None):
    calls = utils.get_intervals_in_region(target_interval, savvy_bed_file,ignore=ignore)
    sample_calls = []
    for call in calls:
        curr_sample = call.data[6].split('.')[0]
        sample_calls.append(call)
    return(sample_calls)

def get_cnvkit_calls(cnvkit_bed_file, target_sample, target_interval):
    calls = utils.get_intervals_in_region(target_interval, cnvkit_bed_file)
    sample_calls = []
    for call in calls:
        cn = int(call.data[1]) 
        if cn == 2:
            continue

        sample = call.data[5]

        if sample != target_sample.split('_')[0]:
            continue

        sample_calls.append(call)
    return(sample_calls)

def man_plot(axs, infile, target):
    cmap= plt.get_cmap('tab10')
    sets = []
    i = 0
    chrm = None
    curr_set = [[],[],[]]
    chrm_labels = [[],[]]
    xs=[]
    means=[]
    stds=[]
    for l in open(infile):
        A = l.rstrip().split('\t')
        if A[0] != target.chrom: continue
        xs.append(int(A[1]))
        means.append(float(A[4]))
        stds.append(float(A[5]))

    axs[0].scatter(xs,means,s=.1,color='#5683d6')
    axs[1].scatter(xs,stds,s=.1,color='#5683d6')

    #axs[0].set_ylabel('Mean', fontsize=6)
    axs[0].set_title('Mean', fontsize=6, x=0.0, fontdict={'horizontalalignment': 'left'})
    #axs[1].set_ylabel('Stdev', fontsize=6)
    axs[1].set_title('Stdev', fontsize=6, x=0.0, fontdict={'horizontalalignment': 'left'})

    axs[1].set_xticks(chrm_labels[0])
    axs[1].set_xticklabels(labels=chrm_labels[1])
    axs[1].axvspan(target.start, target.end,alpha=.5, color='#ff5025')
    axs[0].axvspan(target.start, target.end,alpha=.5, color='#ff5025')

    return axs


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


    samples = utils.get_header(args.scores)[0].split('\t')[4:]
    print(samples)
    print(args.sample)
    sample_i = samples.index(args.sample)
    print(sample_i)
    scores = utils.get_intervals_in_region(target_region,
                                           args.scores)

    exons = utils.get_intervals_in_region(target_region,
                                           args.exons)

    regions_pos = []
    means = []
    stdevs = []
    xs = []
    obs_hets = []
    mean_hets = []
    obs_homs = []
    mean_homs = []
    other = []
    alt_xs = []
    alt_ys = []
    alt_stds = []
    alt_means = []

    non_samp_xs = []    
    non_samp_pos = []
    vcf = VariantFile(args.vcf)

    for exon in scores:
        scores = [float(x) for x in exon.data[1:]]
        regions_pos.append(exon.end)
        means.append(np.mean(scores))
        stdevs.append(np.std(scores))
        xs.append(scores[sample_i])
        non_samp_xs += scores
        non_samp_pos += [exon.end] * len(scores)
        gt_stat =  utils.get_gt_stats(exon, target_sample=args.sample, vcf_handle=vcf)
        if gt_stat:
            obs_homs.append((exon.end,gt_stat[0][1]))
            obs_hets.append((exon.end,gt_stat[1][1]))

            mean_homs.append((exon.end,gt_stat[0][2]))
            mean_hets.append((exon.end,gt_stat[1][2]))
    print(0)
    # find out the number of extra bed files to be included
    num_additional_beds = 0
    if args.additional_beds:
        for line in open(args.additional_beds):
            num_additional_beds += 1
    
    plt.rcParams.update({'font.size': 6})
    axs = []
    default_grid_spec_size = 15
    # add extra axies for extra beds
    extra_bed_increase_factor = num_additional_beds * 2
    grid_spec_size = default_grid_spec_size + extra_bed_increase_factor
    
    fig = plt.figure(constrained_layout=True)
    
    total_num_plots = num_additional_beds + 5 # 5 is the default number of plots
    fig.set_size_inches(6, total_num_plots)

    gs = fig.add_gridspec(grid_spec_size,1, hspace=1.0, wspace=1.0)
    #gs1 = fig.add_gridspec(grid_spec_size, 1, hspace=1.00, wspace=1.00)
    first_chr_plot_index = 3 + num_additional_beds
    num_calls_idx = 1
    gnomad_idx = 2
    last_reg_plt_idx = gnomad_idx + num_additional_beds
    print(1)
    axs.append(fig.add_subplot(gs[:6,:]))
    axs.append(fig.add_subplot(gs[6:8,:]))
    axs.append(fig.add_subplot(gs[8:10,:]))
    if num_additional_beds > 0:
        for i in range(num_additional_beds):
            increase_factor = (i+1) * 2
            axs.append(fig.add_subplot(gs[8+increase_factor:10+increase_factor,:]))
    axs.append(fig.add_subplot(gs[11 + extra_bed_increase_factor:13 + extra_bed_increase_factor,:]))
    axs.append(fig.add_subplot(gs[13 + extra_bed_increase_factor:15 + extra_bed_increase_factor,:]))
    ax_is_discrete = [False] * len(axs) 
    fig.tight_layout()
    for a in axs:
        format_axis(a,args)
    ax = axs[0]
    print(2)
    if args.alt_allele_counts and args.alt_allele_headers:
        tbx = pysam.TabixFile(args.alt_allele_counts)
        for line in open(args.alt_allele_headers):
            t_headers = line.strip().split('\t')
            break
        samp = t_headers.index(args.sample)
        region = args.region.split(':')
        c = region[0]
        s = max(0, int(region[1].split('-')[0]) - args.window)
        e = int(region[1].split('-')[1]) + args.window
        alt_info =  get_region_alt_counts(c,s,e,tbx,samp)
        alt_xs = list(alt_info['start'])
        alt_ys = list(alt_info['zscore'])
        alt_means = list(alt_info['avg'])
        alt_stds = list(alt_info['std'])
        alt_score = list( float(x) for x in alt_info['sample'])
    else:
        alt_score = None
    print(3)
    if alt_score is not None:
        lns10 = ax2.plot(alt_xs,
                   alt_score,
                   '-o',
                   lw = 0,
                   markersize=1,
                   label='Sample Pop. Alt. Count',
                   c='#ff5025')
    else:
        lns10 = None

    ax, lns1, lns2b, lns2, all_calls, savvy_calls, xmin, xmax, target_interval = plot_coverage_stats(ax, args, call_colors, non_samp_pos, non_samp_xs, regions_pos, means, stdevs, xs, chrom)
    print(4)
    if args.gnomad_sv:
        gnomad_sv = pysam.TabixFile(args.gnomad_sv) 
        raw_svs = gnomad_sv.fetch(target_interval.chrom, target_interval.start, target_interval.end)
        sv_xs = []
        sv_ys = []
        for l in raw_svs:
            row = l.strip().split('\t')
            for y in row[37].split(','):
                sv_xs.append(int(row[1]))
                sv_ys.append(float(y))
        axs[gnomad_idx].scatter(sv_xs, sv_ys, s=.1, marker='*')
        #axs[gnomad_idx].set_ylabel('Gnomad SV AF', fontsize=6)
        axs[gnomad_idx].set_title('Gnomad SV AF', fontsize=6, x=0.0, fontdict={'horizontalalignment': 'left'})
    print(5)
    if args.additional_beds:
        for i,line in enumerate(open(args.additional_beds,'r')):
            filename = line.strip()
            basename = filename.split('.')[0]
            extra_bed = pysam.TabixFile(filename)
            bed_entries = extra_bed.fetch(target_interval.chrom, target_interval.start, target_interval.end)
            extra_xs = []
            extra_ys = []
            has_numerics = False
            has_discrete = False
            for entry in bed_entries:
                row = entry.strip().split('\t')
                extra_xs.append(int(row[1]))
                if row[3].isnumeric(): 
                    extra_ys.append(float(row[3]))
                    has_numerics = True
                else:
                    extra_ys.append(1)
                    has_discrete = True
            if has_numerics and has_discrete:
                print('Warning: additional bed', filename, 'has box numeric and non-numeric entries in the 4th column. Treating as all non-numeric')
                extra_ys = [1] * len(extra_xs)
            axs[gnomad_idx + (i+1)].scatter(extra_xs, extra_ys, s=.1, marker='*')
            #axs[gnomad_idx + (i+1)].set_ylabel(basename, fontsize=6)
            axs[gnomad_idx + (i+1)].set_title(basename, fontsize=6, x=0.0, fontdict={'horizontalalignment': 'left'})
            if has_discrete:
                # if they are non-numeric values, hide the y axis ticks
                axs[gnomad_idx + (i+1)].set_yticks([1])
                axs[gnomad_idx + (i+1)].set_yticklabels([' '],)
                axs[gnomad_idx + (i+1)].xaxis.label.set_color('white')
                ax_is_discrete[gnomad_idx + (i+1)] = True
    print(6)
    if args.exons:
        probes = []
        for line in open(args.depth,'r'):
            row = line.strip().split('\t')
            if str(row[0]) != str(target_region.chrom):
                continue
            probes.append(utils.Interval(chrom=row[0],start=int(row[1]),end=int(row[2]),data=None))
        tick_xs = [x.start for x in exons]
        tick_labs = [' '.join(x.data) for x in exons]
        
        axs[last_reg_plt_idx].set_xticks(tick_xs)
        axs[last_reg_plt_idx].set_xticklabels(tick_labs,fontsize= 3,)
        print(7)
        legend_elements = []
        colors = call_colors
        markers = ['o', '^', 's', 'P', 'D', '*']
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
            axs[num_calls_idx].scatter(x_vals, y_vals,s=.1,color=colors[i],marker='*')
            legend_elements.append(Line2D([0], [0],
                              marker='o',
                              color=colors[i],
                              label=get_file_name(f),
                              markerfacecolor=colors[i],
                              markersize=2,
                              linewidth=0))
        axs[first_chr_plot_index].set_xticks([])
        axs[first_chr_plot_index].set_xticklabels([])
        print(8)
        #axs[num_calls_idx].set_ylabel('Num. Calls',fontsize=6)
        axs[num_calls_idx].set_title('Num. Calls', fontsize=6, x=0.0, fontdict={'horizontalalignment': 'left'}) 
        axs[num_calls_idx].tick_params(bottom=False)
        axs[num_calls_idx].legend(handles=legend_elements,frameon=False,prop={'size': 5},bbox_to_anchor=(1.10,1), borderaxespad=0)
        if max(density_of_calls_at_each_probe) < 3 and False:
            print('Changed!!!')
            ticks = list(range(max(density_of_calls_at_each_probe)+1))
            axs[num_calls_idx].set_yticks(ticks)
            axs[num_calls_idx].set_yticklabels(ticks)
        else:
            print('no change!!')
        print(9)
        if args.max_num_calls:
            # process this file
            mnc = process_max_num_calls(args.max_num_calls)
            axs[num_calls_idx].set_ylim((axs[num_calls_idx].get_ylim()[0],mnc))
    print(10)
    lns_to_include = [lns1, lns2, lns10, lns2b]
    lns = None
    for x in lns_to_include:
        if x is not None:
            if lns is None:
                lns = x
            else:
                lns += x
    
    man_plot([axs[first_chr_plot_index],axs[first_chr_plot_index+1]], args.depth, target_region)
    labs = [l.get_label() for l in lns]
    ax.set_ylim([-7, 7])
    leg = ax.legend(lns,
                    labs,
                    frameon=False,
                    prop={'size': 5},bbox_to_anchor=(1.10,1), borderaxespad=0)
    print(11)
    for i in  range(0,len(axs)):
        axs[i].tick_params(axis='both', which='major', labelsize=4,length=2,width=.5)
        axs[i].tick_params(axis='both', which='major', labelsize=4,length=2,width=.5)

    axs[last_reg_plt_idx].tick_params(axis='both', which='major', labelsize=4,length=2,width=.5)
    axs[last_reg_plt_idx].tick_params(axis='both', which='minor', labelsize=4,length=2,width=.5)

    axs[first_chr_plot_index].set_xticks([])
   
    # remove labels on all regional plots
    for i in range(last_reg_plt_idx):
        axs[i].set_xticks([])
        axs[i].set_xticklabels([])

    # add exon labels to the last regional plot
    axs[last_reg_plt_idx].tick_params(axis='x', which='major', labelsize=2, rotation=90)
    axs[last_reg_plt_idx].tick_params(axis='x', which='minor', labelsize=2, rotation=90)

    for i in range(len(axs)):
        axs[i].spines['top'].set_visible(False)
        axs[i].spines['right'].set_visible(False)
        axs[i].spines['bottom'].set_visible(False)
        axs[i].spines['left'].set_visible(False)
    axs[-1].spines['left'].set_visible(False)
    axs[-1].set_xlabel('Chromosome ' + str(target_region.chrom), fontsize=6)
    axs[last_reg_plt_idx].spines['bottom'].set_visible(True)

    regional_min = min(ax.get_xlim()[0],ax.get_xlim()[0])
    regional_max = max(ax.get_xlim()[1],ax.get_xlim()[1])

    mins = []
    maxs = []
    print(12)
    # ensure the chromosomal plots have the same x limits
    for i in range(first_chr_plot_index, len(axs)):
        mins.append(axs[i].get_xlim()[0])
        maxs.append(axs[i].get_xlim()[1])
    for i in range(first_chr_plot_index, len(axs)):
        axs[i].set_xlim((min(mins), max(maxs)))

    # make sure that all of the regional plots have the same xlimits
    for i in range(last_reg_plt_idx+1):
        axs[i].set_xlim((regional_min,regional_max))
    print(13)
    # rasterize all the points
    for i in range(len(axs)):
        axs[i].set_rasterized(True)
        if i > 0 and i < first_chr_plot_index:
            if not ax_is_discrete[i]:
                axs[i].grid(axis='y',alpha=.2)
    print(14)
    plt.gca().spines['right'].set_color('none')
    plt.gca().spines['top'].set_color('none')
    plt.gca().spines['left'].set_color('none')
    plt.savefig(args.outfile,bbox_inches='tight',dpi=600)

if __name__ == '__main__': main()
