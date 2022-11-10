#!/usr/bin/env python
import argparse
import sys
import os

def get_args():
    help_message ="""
        SeeNV usage:

        One of the following options is required

        --buildref, -b      flag to use function to build a reference panel database
        --plotsamples, -p   flag to plot samples

        Parameters to accompany --plotsamples, -p:
            -i INPUT_SAMPLES   (required) samples list
            -s SITES           (required) genomic sites bed file
            -c ALL_CALLS       (required) calls file. Each line should be a path to a set of calls in bed format
            -o OUTPUT          (required) output directory, where to save the plots
            -r REF_DB          (required) path to reference db created by the --buildref function of SeeNV
            -g GENES_FILE      (required) a bed format genes file
            -a GNOMAD          (required) the gnomad sv sites file with allele frequency information
            -t THREADS         (optional) number of threads to use, default 1 (you really want to use more than 1)

        Parameters to accompany --buildref, -b
            -i INPUT_SAMPLES  (required) samples list
            -s SITES          (required) genomic sites bed file
            -c ALL_CALLS      (required) calls file. Each line should be a path to a set of calls in bed format
            -o OUTPUT         (required) output directory, reference panel database
            -t THREADS        (optional) number of threads to use, default 1 (you really want to use more than 1)
        """
    
    if '-h' in sys.argv or '--help' in sys.argv:
        print(help_message)
    param_sum = sum(['--buildref' in sys.argv or '-b' in sys.argv, '--plotsamples' in sys.argv or '-p' in sys.argv])
    if param_sum > 1 or param_sum == 0:
        print('SeeNV Error: you must use --buildref (-b) xor --plotsamples (-p)')
        print(help_message)
        exit(1)
    run_type = None
    try:
        if '--buildref' in sys.argv or '-b' in sys.argv:
            args = get_build_args()
            run_type = 'buildref'
        elif '--plotsamples' in sys.argv or '-p' in sys.argv:
            args = get_plot_args()
            run_type = 'plotsamples'
    except:
        print(help_message)
        exit(1)
    return run_type, args

def get_plot_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plotsamples',
    '-p',
    action='store_true',
    dest='plotsamples',
    default=False,
    help='flag to plot samples')
    parser.add_argument('-i',dest='input_samples',help='samples list',required=True)
    parser.add_argument('-s',dest='sites',help='genomic sites bed file',required=True)
    parser.add_argument('-c',dest='all_calls',help='calls file. Each line should be a path to a set of calls in bed format',required=True)
    parser.add_argument('-o',dest='output',help='output directory, where to save the plots',required=True)
    parser.add_argument('-r',dest='ref_db',help='path to reference db created by the --buildref function of SeeNV',required=True)
    parser.add_argument('-g',dest='genes_file',help='a bed format genes file',required=True)
    parser.add_argument('-a',dest='gnomad',help='the gnomad sv sites file with allele frequency information',required=True)
    parser.add_argument('-t',dest='threads',help='number of threads to use, default 1 (you really want to use more than 1)',default="1",required=False)
    return parser.parse_args()

def get_build_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--buildref',
                        '-b',
                        action='store_true',
                        dest='buildref',
                        default=False,
                        help='flag to use function to build a reference panel database')
    parser.add_argument('-i',dest='input_samples',help='samples list',required=True)
    parser.add_argument('-s',dest='sites',help='genomic sites bed file',required=True)
    parser.add_argument('-c',dest='all_calls',help='calls file. Each line should be a path to a set of calls in bed format',required=True)
    parser.add_argument('-o',dest='output',help='output directory, reference panel database',required=True)
    parser.add_argument('-t',dest='threads',help='number of threads to use, default 1 (you really want to use more than 1)',default="1",required=False)
    return parser.parse_args()

run_type, args = get_args()
if run_type == 'plotsamples':
    command="""
    conda_loc=$(which python)
    conda_bin=$(dirname $conda_loc)
    inputsamples="{samples}"
    probes="{probes}"
    calls="{calls}"
    genesfile="{genes_file}"
    gnomadfile="{af}"
    outputdir="{out}"
    refpanel="{ref_db}"
    threads="{threads}"

    mkdir -p workproband
    cp $inputsamples workproband/proband.samples
    cat $inputsamples | bin/gargs --sep="\t" "ln -f -s {{3}} workproband/{{0}}.bam"
    cat $inputsamples | bin/gargs --sep="\t" "ln -f -s {{4}} workproband/{{0}}.bai"
    cat $inputsamples | bin/gargs --sep="\t" "ln -f -s {{5}} workproband/{{0}}.vcf.gz"
    cat $inputsamples | bin/gargs --sep="\t" "ln -f -s {{6}} workproband/{{0}}.vcf.gz.tbi"

    mkdir -p workproband/Probes
    ln -s "$(cd "$(dirname "$probes")"; pwd)/$(basename "$probes")" workproband/Probes/probes.original.bed
    mkdir -p workproband/Calls
    ln -s "$(cd "$(dirname "$calls")"; pwd)/$(basename "$calls")" workproband/Calls/all.calls.original.txt
    ln -s "$(cd "$(dirname "$genesfile")"; pwd)/$(basename "$genesfile")" workproband/genes.bed.gz
    ln -s "$(cd "$(dirname "$gnomadfile")"; pwd)/$(basename "$gnomadfile")" workproband/gnomad_sv.bed.gz
    ln -s "$(cd "$(dirname "$gnomadfile")"; pwd)/$(basename "$gnomadfile")".tbi workproband/gnomad_sv.bed.gz.tbi

    # put the name of the output dirrectory into a textfile for Snakemake to read in
    echo $outputdir > workproband/outputdir.txt
    # put the path to the reference panel into a textfie for Snakemake to read in
    echo $refpanel > workproband/reference_panel.txt
    
    # start the snakemake pipeline now they input files are in there proper locations
    snakemake -c $threads -s $conda_bin/run_proband.snake --configfile $conda_bin/proband_config.json --printshellcmds --rerun-incomplete --restart-times 3
    """
    command = command.format(samples=args.input_samples,
                    probes=args.sites,
                    calls=args.all_calls,
                    out=args.output,
                    ref_db=args.ref_db,
                    genes_file=args.genes_file,
                    af=args.gnomad,
                    threads=args.threads)
elif run_type == 'buildref':
    print('Building')
    command="""
    conda_loc=$(which python)
    conda_bin=$(dirname $conda_loc)

    inputsamples="{samples}"
    probes="{probes}"
    calls="{calls}"
    outputdir="{out}"
    threads="{threads}"

    # create the workpanel dir, dump sym linked files in to a location they can be easily accessed by SnakeMake 
    # will forably overwrite input files of the same name. All samples named MUST be unique.
    mkdir -p workpanel
    cp $inputsamples workpanel/panel.samples
    cat $inputsamples | bin/gargs --sep="\t" "ln -f -s {{3}} workpanel/{{0}}.bam"
    cat $inputsamples | bin/gargs --sep="\t" "ln -f -s {{4}} workpanel/{{0}}.bai"
    cat $inputsamples | bin/gargs --sep="\t" "ln -f -s {{5}} workpanel/{{0}}.vcf.gz"
    cat $inputsamples | bin/gargs --sep="\t" "ln -f -s {{6}} workpanel/{{0}}.vcf.gz.tbi"
 
    mkdir -p workpanel/Probes
    ln -f -s "$(cd "$(dirname "$probes")"; pwd)/$(basename "$probes")" workpanel/Probes/probes.original.bed
    mkdir -p workpanel/Calls
    ln -f -s "$(cd "$(dirname "$calls")"; pwd)/$(basename "$calls")" workpanel/Calls/all.calls.original.txt
    
    # put the name of the output dirrectory into a textfile for Snakemake to read in
    echo $outputdir > workpanel/outputdir.txt

    # start the snakemake pipeline now they input files are in there proper locations
    snakemake -c $threads -s $conda_bin/build_panel.snake --configfile $conda_bin/panel_config.json --restart-times 3
    """
    command = command.format(samples=args.input_samples,probes=args.sites,calls=args.all_calls,out=args.output,threads=args.threads)
else:
    print('Internal Error, unexpected run_type')
    exit(1)

os.system(command)
