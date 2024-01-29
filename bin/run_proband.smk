import pandas as pd
import os

tempdir = 'workproband'
if 'workdir' in config:
	tempdir = config['workdir']
if tempdir[-1] != '/':
	tempdir += '/'

# ensure probes, all_calls gnomad_sv_bed_gz gnomad_sv_bed_gz_tbi vardb_bed_gz repeat_masker_gz all have tempdir in the front of them
for key in ['probes','all_calls','gnomad_sv_bed_gz','gnomad_sv_bed_gz_tbi','vardb_bed_gz','repeat_masker_gz']:
	if  tempdir not in config[key]:
		config[key] = tempdir + config[key]

# Load input files and options
samples_df = pd.read_csv(tempdir + 'proband.samples',sep='\t',header=None)
samples_df.columns = ['p_id','s_id','sex','bam','bai']

SAMPLES = samples_df['s_id']
BAM = samples_df["bam"]

for line in open(tempdir + 'outputdir.txt','r'):
	config['outputdir'] = line.strip()
	break

for line in open(tempdir + 'reference_panel.txt','r'):
	config['reference_panel'] = line.strip()
	break

# get full absolute path to the reference panel
config['reference_panel'] = os.path.abspath(config['reference_panel'])
if config['reference_panel'][-1] != '/':
	config['reference_panel'] += '/'

# get all reference panel RPM files
REF_RPMs = []
for f in os.listdir(config['reference_panel'] + 'RPM/'):
	REF_RPMs.append(config['reference_panel'] + 'RPM/' + f)

ref_samples_df = pd.read_csv(config['reference_panel'] + 'panel.samples',sep='\t',header=None)
ref_samples_df.columns = ['p_id','s_id','sex','bam','bai']
REF_SAMPLES = list(ref_samples_df['s_id'])
	
num_calls = 0
for line in open(config['all_calls'],'r'):
	if sum( s in line for s in SAMPLES) > 0:
		num_calls += 1
NUM_CALLS = [str(x) for x in range(num_calls)]

if len(NUM_CALLS) == 0:
	print('SeeNV Error: No calls found for the samples in the all_calls file')
	exit(1)

def get_dir_name(filepath):
	"""
	Get the directory name from a file path
	"""
	name = os.path.dirname(filepath)
	if name[-1] != '/':
		name += '/'
	return name

rule all:
	input:
		expand(tempdir + "Mosdepth/{sample}.per-base.bed.gz", sample=SAMPLES),
		tempdir + "Probes/probes.sorted.bed.gz",
		tempdir + "Probes/probes.sorted.bed.gz.tbi",
		expand(tempdir + "ReadCounts/{sample}.num_reads.txt", sample=SAMPLES),
		tempdir + "TotalReads/total_read.txt",
		expand(tempdir + "RPM/{sample}.probe.rpm_rate.bed.gz", sample=SAMPLES),
		expand(tempdir + "RPM/{sample}.probe.rpm_rate.bed.gz.tbi", sample=SAMPLES),
		tempdir + "Params/call_plotting_params.txt",
		#tempdir + "Exons/labeled_exons.bed.gz",
		#tempdir + "Exons/labeled_exons.bed.gz.tbi",
		expand(tempdir + "ProbeCoverage/{sample}_probe.cover.mean.stdev.bed.gz", sample=SAMPLES),
		expand(tempdir + "ProbeCoverage/{sample}_probe.cover.mean.stdev.bed.gz.tbi", sample=SAMPLES),
		expand(tempdir + 'ProbeCoverageRefPanel/{ref_sample}.probe.rpm_rate.bed.gz', ref_sample=REF_SAMPLES),
		expand(tempdir + 'ProbeCoverageRefPanel/{ref_sample}.probe.rpm_rate.bed.gz.tbi', ref_sample=REF_SAMPLES),
		expand(tempdir + "AdjZscore/{sample}.adj_z.bed.gz", sample=SAMPLES),
		expand(tempdir + "AdjZscore/{sample}.adj_z.bed.gz.tbi", sample=SAMPLES),
		tempdir + "MergedAdjZscore/adj_scores.bed.gz",
		tempdir + "MergedAdjZscore/adj_scores.bed.gz.tbi",
		tempdir + "SplitParams/.split_param_file.ready",
		expand(tempdir + "PlotsComplete/{num}.done",num=NUM_CALLS),
		expand(tempdir + "AdjZscoreRefPanel/{ref_sample}.adj_z.bed.gz", ref_sample=REF_SAMPLES),
		expand(tempdir + "AdjZscoreRefPanel/{ref_sample}.adj_z.bed.gz.tbi", ref_sample=REF_SAMPLES),
		expand(tempdir + "AlleleCount/{sample}.allele_count.tsv", sample=SAMPLES),
		tempdir + "all.done"


rule mosdepth:
	input:
		tempdir + "{sample}.bam"
	output:
		tempdir + "Mosdepth/{sample}.per-base.bed.gz",
		tempdir + "Mosdepth/{sample}.per-base.bed.gz.tbi"
	log:
		tempdir + "logs/mosdepth.{sample}.log"
	params:
		tempdir = tempdir,
		outputdirname = tempdir + "Mosdepth/"
	shell:
		"""
		mkdir -p {params.outputdirname}
		mosdepth {params.outputdirname}{wildcards.sample} {input}
		tabix -f -p bed {params.outputdirname}{wildcards.sample}.per-base.bed.gz
		"""

rule count_reads:
	input:
			tempdir + "{sample}.bam"
	output:
			tempdir + "ReadCounts/{sample}.num_reads.txt"
	log:
			tempdir + "logs/count_reads.{sample}.log"
	params:
		outputdirname = tempdir + "ReadCounts/"
	shell:
		"""
		mkdir -p {params.outputdirname}
		samtools view -c -F 260 {input} > {output}
		"""

rule get_total_read_counts:
#	"""
#	Add all the individual read counts to one file	
#	"""
	input:
		expand(tempdir + "ReadCounts/{sample}.num_reads.txt", sample=SAMPLES)
	output:
		tempdir + "TotalReads/total_read.txt"
	params:
		outputdirname = tempdir + "TotalReads/"
	shell:
		"""
		mkdir -p {params.outputdirname}
		cat {input} > {params.outputdirname}total_read.txt
		"""

rule gzip_probes:
#	"""
#	Take the probes files, bgzip and tabix it
#	"""
	input:
		probes = config['probes']
	output:
		tempdir + "Probes/probes.sorted.bed.gz",
		tempdir + "Probes/probes.sorted.bed.gz.tbi"
	params:
		outputdirname = tempdir + "Probes/"
	shell:
		"""
		mkdir -p {params.outputdirname}
		bedtools sort -i {input.probes} > {params.outputdirname}probes.sorted.bed
		bgzip {params.outputdirname}probes.sorted.bed
		tabix {params.outputdirname}probes.sorted.bed.gz -p bed
		"""

rule get_probe_reads_per_million:
#	"""
#	For each sample calculate the number of reads per million (RPM) at each probe. Return the files bgzip-ed and tabix-ed
#	"""
	input:
		per_base_bed_gz=tempdir + "Mosdepth/{sample}.per-base.bed.gz",
		total_reads=tempdir + "TotalReads/total_read.txt",
		probes_gz=tempdir + "Probes/probes.sorted.bed.gz"
	output:
		tempdir + "RPM/{sample}.probe.rpm_rate.bed.gz",
		tempdir + "RPM/{sample}.probe.rpm_rate.bed.gz.tbi"
	params:
		outputdirname = tempdir + "RPM/"
	shell:
		"""
		mkdir -p {params.outputdirname}
		get_rpm_rates.py \
		-m {input.per_base_bed_gz} \
		--regions_file {input.probes_gz} \
		--num_reads {input.total_reads} \
		| bgzip -c > {params.outputdirname}{wildcards.sample}.probe.rpm_rate.bed.gz
		tabix -p bed {params.outputdirname}{wildcards.sample}.probe.rpm_rate.bed.gz
		"""

rule create_plotting_params:
#	"""
#	Create a file will add the CNVs and their preconficgures parameters needed for plotting
#	"""
	input:
		all_calls=config['all_calls']

	output:
		tempdir + "Params/call_plotting_params.txt"
	log:
		tempdir + "logs/plotting_params.log"
	params:
		outputdirname = tempdir + "Params/",
		samples=expand("{sample}",sample=SAMPLES)
	shell:
		"""
		mkdir -p {params.outputdirname}
		cat {input.all_calls} | awk '{{ print $1":"$2"-"$3"\t"$4"\t"$5"\t.\t"$1"."$2"-"$3"\t"$4}}' >> {params.outputdirname}tmp_call_plotting_params.txt
		
		for i in {{1..23}}; do
			grep "^$i" {params.outputdirname}tmp_call_plotting_params.txt >> {params.outputdirname}call_plotting_params.txt.tmp || echo "Not found" 2>> {log}
		done

		touch {params.outputdirname}call_plotting_params.txt
		# remove samples from the calls that are not in the sample list
		array=( {params.samples} )
		for j in "${{array[@]}}"
		do
			cat {params.outputdirname}call_plotting_params.txt.tmp | {{ grep $j || :; }} >> {params.outputdirname}call_plotting_params.txt
		done
		"""

checkpoint split_param_file:
#	"""
#	Take the one large param files and split it up into files with a single line in it
#	"""
	input:
		param_file=tempdir + "Params/call_plotting_params.txt",
	output:
		tempdir + "SplitParams/.split_param_file.ready"
	log:
		tempdir + "logs/split_param_file.log"
	params:
		outputdirname = tempdir + "SplitParams/",
		samples=expand("{sample}",sample=SAMPLES)
	shell:
		"""
		mkdir -p {params.outputdirname}
		i=0
		while read p; do
			echo $p > {params.outputdirname}${{i}}.txt
			i=$((i+1))
		done <{input.param_file}
		touch {params.outputdirname}.split_param_file.ready
		"""

class Checkpoint_MakePattern_prime:
	"""
	Class that makes dealing with checkpoints easy. Take from:
	http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
	"""
	def __init__(self, pattern):
			self.pattern = pattern

	def __call__(self, w):
		global checkpoints

		checkpoints.split_param_file.get(**w)

		names = glob_wildcards(tempdir + 'SplitParams/{rs}.txt')[0]

		return self.pattern

rule label_exons:
	"""
	Take in the genes files and add create named labels for the exons in all genes.
	Return it bgzip-ed and tabix-ed
	"""
	input:
		tempdir + "genes.bed.gz"
	output:
		tempdir + "Exons/labeled_exons.bed.gz",
		tempdir + "Exons/labeled_exons.bed.gz.tbi"
	params:
		outputdirname = tempdir + "Exons/"
	shell:
		"""
		mkdir -p {params.outputdirname}
		zcat {input} | label_exon_number.py > {params.outputdirname}labeled_exons.bed
		bgzip {params.outputdirname}labeled_exons.bed
		tabix -p bed {params.outputdirname}labeled_exons.bed.gz
		"""

rule link_ref_panel_rpm:
	"""
	Make a sym-link of the reference panels RPM files for later use with each proband
	"""
	input:
		bedgz=config['reference_panel'] + "RPM/{ref_sample}.probe.rpm_rate.bed.gz",
		tbi=config['reference_panel'] + "RPM/{ref_sample}.probe.rpm_rate.bed.gz.tbi",
		adj_bed=config['reference_panel'] + "AdjZscore/{ref_sample}.adj_z.bed.gz",
		adj_tbi=config['reference_panel'] + "AdjZscore/{ref_sample}.adj_z.bed.gz.tbi"	
	output:
		tempdir + 'ProbeCoverageRefPanel/{ref_sample}.probe.rpm_rate.bed.gz',
		tempdir + 'ProbeCoverageRefPanel/{ref_sample}.probe.rpm_rate.bed.gz.tbi',
		tempdir + 'AdjZscoreRefPanel/{ref_sample}.adj_z.bed.gz',
		tempdir + 'AdjZscoreRefPanel/{ref_sample}.adj_z.bed.gz.tbi'

	params:
		outputdirname = tempdir + 'ProbeCoverageRefPanel/',
		outputdirname2 = tempdir + 'AdjZscoreRefPanel/'
	shell:
		"""
		mkdir -p {params.outputdirname}
		mkdir -p {params.outputdirname2}

		ln -f -s {input.bedgz} {params.outputdirname}
		ln -f -s {input.tbi} {params.outputdirname}
		ln -f -s {input.adj_bed} {params.outputdirname2}
		ln -f -s {input.adj_tbi} {params.outputdirname2}
		"""


rule get_probe_cover_mean_std:
	"""
	Calculate various stats about RPM coverage by comparing one proband to the reference panel
	"""
	input:
		ref_rpm=expand(tempdir + 'ProbeCoverageRefPanel/{ref_sample}.probe.rpm_rate.bed.gz', ref_sample=REF_SAMPLES),
		ref_rpm_tbi=expand(tempdir + 'ProbeCoverageRefPanel/{ref_sample}.probe.rpm_rate.bed.gz.tbi', ref_sample=REF_SAMPLES),
		proband_rpm=tempdir + "RPM/{sample}.probe.rpm_rate.bed.gz"	
	output:
		bedgz=tempdir + "ProbeCoverage/{sample}_probe.cover.mean.stdev.bed.gz",
		tbi=tempdir + "ProbeCoverage/{sample}_probe.cover.mean.stdev.bed.gz.tbi",
		bed=tempdir + "ProbeCoverage/{sample}_probe.cover.mean.stdev.bed"
	threads: 2
	params:
		outputdirname = tempdir + "ProbeCoverage/",
	shell:
		"""
		mkdir -p {params.outputdirname}
		get_regions_zscores.py -r {params.outputdirname} -s {input.proband_rpm} | bgzip -c > {output.bedgz}
		cp {output.bedgz} {params.outputdirname}{wildcards.sample}_probe.cover.mean.stdev.copy.bed.gz
		gunzip {params.outputdirname}{wildcards.sample}_probe.cover.mean.stdev.copy.bed.gz
		mv {params.outputdirname}{wildcards.sample}_probe.cover.mean.stdev.copy.bed {output.bed}
		tabix -p bed {output.bedgz}
		"""


rule get_adj_zscore:
#	"""
#	Calculate the z-score of the RPM data compared to the std and mean for each proband 
#	"""
	input:
		proband_rpm=tempdir + "RPM/{sample}.probe.rpm_rate.bed.gz",
		proband_coverage=tempdir + "ProbeCoverage/{sample}_probe.cover.mean.stdev.bed.gz",
		proband_coverage_tbi=tempdir + "ProbeCoverage/{sample}_probe.cover.mean.stdev.bed.gz.tbi"
	output:
		bedgz=tempdir + "AdjZscore/{sample}.adj_z.bed.gz",
		tbi=tempdir + "AdjZscore/{sample}.adj_z.bed.gz.tbi"
	params:
		outputdirname = tempdir + "AdjZscore/",
	shell:
		"""
		mkdir -p {params.outputdirname}
		get_coverage_zscores.py \
			-r {input.proband_rpm} \
			-s {input.proband_coverage} \
			| bgzip -c > {output.bedgz}
		tabix -p bed {output.bedgz}	
		"""

rule merge_adj_scores:
	"""
	Take the adj_z_scores for every sample and put them into a single file
	"""
	input:
		bedgz=expand(tempdir + "AdjZscore/{sample}.adj_z.bed.gz", sample=SAMPLES),
		tbi=expand(tempdir + "AdjZscore/{sample}.adj_z.bed.gz.tbi", sample=SAMPLES),
		ref_bed=expand(tempdir + "AdjZscoreRefPanel/{ref_sample}.adj_z.bed.gz", ref_sample=REF_SAMPLES),
		ref_tbi=expand(tempdir + "AdjZscoreRefPanel/{ref_sample}.adj_z.bed.gz.tbi", ref_sample=REF_SAMPLES)
	output:
		final_bed_gz=tempdir + "MergedAdjZscore/adj_scores.bed.gz",
		final_bed_gz_tbi=tempdir + "MergedAdjZscore/adj_scores.bed.gz.tbi"
	params:
		outputdirname = tempdir + "MergedAdjZscore/",
		tempdir = tempdir
	shell:
		"""
		mkdir -p {params.outputdirname}
		mkdir -p {params.tempdir}AllAdjZScore
		cp {params.tempdir}AdjZscore/* {params.tempdir}AllAdjZScore
		cp {params.tempdir}AdjZscoreRefPanel/* {params.tempdir}AllAdjZScore
		merge_samples_adj_scores.py -r {params.tempdir}AllAdjZScore | bgzip -c > {output.final_bed_gz}
		tabix -p bed {output.final_bed_gz}
		"""

rule prep_all_calls:
	input:
		all_calls=config['all_calls']
	output:
		tempdir + "FindMaxTMP/.ready"
	log:
		tempdir + "logs/prep_all_calls.log"
	params:
		outputdirname = tempdir + "FindMaxTMP/",
		tempdir = tempdir,
		samples=expand("{sample}",sample=SAMPLES)
	shell:
		"""
		mkdir -p {params.outputdirname}
		rm -f {params.outputdirname}*
		
		name=$(basename {input.all_calls})
		array=( {params.samples} )
		for i in "${{array[@]}}"
		do
			cat {input.all_calls} | {{ grep $i || :; }} >> {params.outputdirname}$name.filtered 2>> {log}
		done

		bedtools sort -i {params.outputdirname}$name.filtered > {params.outputdirname}$name.sorted 2>> {log}
		bgzip -f {params.outputdirname}$name.sorted 2>> {log}
		tabix -f -p bed {params.outputdirname}$name.sorted.gz 2>> {log}
		echo "{params.outputdirname}${{name}}.sorted.gz" >> {params.outputdirname}multiple_savvy_calls.txt 2>> {log}
		
		touch {params.outputdirname}.ready
		"""

rule plotter:
	input:
		probes=tempdir + "Probes/probes.sorted.bed.gz",
		probes_tbi=tempdir + "Probes/probes.sorted.bed.gz.tbi",
		all_calls=tempdir + "FindMaxTMP/.ready",
		param_file=Checkpoint_MakePattern_prime(tempdir + 'SplitParams/{num}.txt'),
		marker_file=tempdir + "SplitParams/.split_param_file.ready",
		coverage=expand(tempdir + "ProbeCoverage/{sample}_probe.cover.mean.stdev.bed.gz", sample=SAMPLES),
		coverage_tbi=expand(tempdir + "ProbeCoverage/{sample}_probe.cover.mean.stdev.bed.gz.tbi", sample=SAMPLES),
		adj_scores=tempdir + "MergedAdjZscore/adj_scores.bed.gz",
		adj_scores_tbi=tempdir + "MergedAdjZscore/adj_scores.bed.gz.tbi",
		gnomad_sv=config['gnomad_sv_bed_gz'],
		vardb=config['vardb_bed_gz'],
		repeat_masker=config['repeat_masker_gz']

	output:
		tempdir + "PlotsComplete/{num}.done"
	log:
		tempdir + "logs/plotting_{num}.log"
	params:
		num_calls=len(NUM_CALLS),
		site_quality=config['site_quality'],
		outputdirname = tempdir + "PlotsComplete/",
		tempdir = tempdir
	shell:
		"""
		mkdir -p {params.outputdirname}
		mkdir -p {params.tempdir}Plots/
		inputs=""
		empty=""
		while read p; do
			if [ "$inputs" == "$empty" ]; then
				inputs="$p"
			else
				inputs="$inputs,$p"
			fi
		done < {params.tempdir}FindMaxTMP/multiple_savvy_calls.txt
		echo {params.num_calls} >> {log}
		echo {wildcards.num} >> {log}
		mkdir -p {params.tempdir}Plots/ 2>> {log}
		mkdir -p {params.tempdir}PlotsComplete/ 2>> {log}
		string=$(head -1 {input.param_file}) 2>> {log}

		region=$(cut -d" " -f1 <<< $string) 2>> {log}
		svtype=$(cut -d" " -f2 <<< $string) 2>> {log}
		samplename=$(cut -d" " -f3 <<< $string) 2>> {log}
		sample=${{samplename%.*}} 2>> {log}
		callsfile=$(cut -d" " -f4 <<< $string) 2>> {log}
		regionclean=$(cut -d" " -f5 <<< $string) 2>> {log}
		svtype2=$(cut -d" " -f6 <<< $string) 2>> {log}

		plotter.py \
			-i $inputs \
			-s ${{samplename%.*}} \
			--output {params.tempdir}Plots/"${{sample}}.${{region}}.${{svtype}}".png \
			--coverage_scores {input.adj_scores} \
			--sites {input.probes} \
			--window 50000 \
			--region ${{region}} \
			--title "${{sample}} ${{region}} ${{svtype}}" \
			--depth {params.tempdir}ProbeCoverage/${{sample}}_probe.cover.mean.stdev.bed \
			--gnomad {input.gnomad_sv} \
			--vardb {input.vardb} \
			--repeatmasker {input.repeat_masker} \
			--site_quality {params.site_quality}
		touch {output}
		"""

rule collect_allele_counts:
	input:
		bam=tempdir + "{sample}.bam",
		probes=config['probes']
	output:
		tempdir + "AlleleCount/{sample}.allele_count.tsv"
	resources:
		mem_mb=20000
	params:
		outputdirname = tempdir + "AlleleCount/",
		tempdir = tempdir
	shell:
		"""
		mkdir -p {params.outputdirname}
		touch {output}
		"""

rule collect_plots:
	input:
		expand(tempdir + "PlotsComplete/{num}.done",num=NUM_CALLS)
	output:
		indicator=tempdir + "all.done"
	params:
		outputdir=config['outputdir'],
		tempdir = tempdir
	shell:
		"""
		mkdir -p {params.outputdir}
		cp {params.tempdir}Plots/* {params.outputdir}
		touch {output.indicator}
		"""