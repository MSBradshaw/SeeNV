import argparse
import sys
import os

def get_args():
    help_message ="""
        CNViz usage:

        One of the following options is required

        --buildref, -b      flag to use function to build a reference panel database
        --plotsamples, -p   flag to plot samples

        Parameters to accompany --plotsamples, -p:
            -i INPUT_SAMPLES   (required) samples list
            -s SITES           (required) genomic sites bed file
            -c ALL_CALLS       (required) calls file. Each line should be a path to a set of calls in bed format
            -o OUTPUT          (required) output directory, where to save the plots
            -r REF_DB          (required) path to reference db created by the --buildref function of CNViz
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
        print('CNViz Error: you must use --buildref (-b) xor --plotsamples (-p)')
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
    parser.add_argument('-r',dest='ref_db',help='path to reference db created by the --buildref function of CNViz',required=True)
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
    return parser.parse_args()

run_type, args = get_args()
if run_type == 'plotsamples':
    print('Plotting')
    command='./start_proband.sh -i {samples}  -p {probes} -c {calls} -o {out} -r {ref_db} -g {genes_file} -a {af} -t {threads}'
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
else:
    print('Internal Error, unexpected run_type')
    exit(1)

os.system(command)
