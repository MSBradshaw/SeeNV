from statsmodels.distributions.empirical_distribution import ECDF
import pysam
import sys
import argparse
import pickle

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--all_rpms',
                        dest='all_rpms',
                        nargs="+",
                        help='space seporated list of RPM filepaths to create ECDFs for each probe with')
    parser.add_argument('--probes',dest='probes',help='path to bed formated list of probes')
    parser.add_argument('--compare_rpms',nargs="+",dest='compare_rpms',help='space seporated list of RPM filepaths to be compared to the ECDF')
    args = parser.parse_args()
    return args
args = get_args()

# reference RPM files
# a single RPM file to compare against
# load all of the RPM files
if os.path.exists('probe_rpms.pickle') and os.path.exists('probe_ecdfs.pickle'):
    probe_rpms = pickle.load(open('probe_rpms.pickle','rb'))
    probe_ecdfs = pickle.load(open('probe_ecdfs.pickle','rb'))
else:
    tbx_objects = []
    for f in args.all_rpms:
        tbx = pysam.TabixFile(f)
        tbx_objects.append(tbx)
    print('Number of tabix objects: ', str(len(tbx_objects)))

# for each probe

    probe_ecdfs = {}
    total_errors = 0
    probe_rpms = {}
    for line in open(args.probes, 'r'):
        row = line.strip().split('\t')
        chrom = row[0]
        start = int(row[1])
        end = int(row[2])
        name = row[3]
        # create an ECDF of each probe
        this_probes_rpms = []
        error_count = 0
        for tbx in tbx_objects:
            x = tbx.fetch(chrom,start,end)
            # for each line in the RPM file that match the query (should be exactly 1)
            match_count = 0
            for l in x:
                r = l.strip().split('\t')
                if r[0] == chrom and start == int(r[1]) and end == int(r[2]):
                    match_count += 1
                    this_probes_rpms.append(float(r[4]))
            if match_count != 1:
                #print('FUCK! There is not exactly 1 match for the probe query in this sample')
                error_count += 1
        probe_ecdfs[name] = ECDF(this_probes_rpms)
        probe_rpms[name] = this_probes_rpms
        print('Error count: ',str(error_count))
        total_errors += error_count
    print('Final Error Count',str(total_errors))
        
    pickle.dump(probe_rpms, open('probe_rpms.pickle','wb'))
    pickle.dump(probe_ecdfs, open('probe_ecdfs.pickle','wb'))

for f in args.compare_rpms:
    for line in open(f,'r'):
        row = line.strip().split()




