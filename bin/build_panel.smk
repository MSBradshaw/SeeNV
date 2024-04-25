import pandas as pd

tempdir = 'workpanel'
if 'workdir' in config:
	tempdir = config['workdir']
if tempdir[-1] != '/':
	tempdir += '/'

# for config items probes, all_calls, samples, probes ensure they all start with tempdir
for key in ['probes', 'all_calls', 'samples']:
	if key in config:
		if config[key].startswith(tempdir) == False:
			config[key] = tempdir + config[key]

SAMPLES_FILE = config['samples']

# Load input files and options
samples_df = pd.read_csv(SAMPLES_FILE,sep='\t',header=None)
samples_df.columns = ['p_id','s_id','sex','bam','bai']

PROBESFILE = config['probes']

SAMPLES = samples_df['s_id']
BAM = samples_df["bam"]

if config['outputdir'][-1] != '/':
	config['outputdir'] += '/'

print(config)

rule all:
	input:
		expand(tempdir + "Mosdepth/{sample}.per-base.bed.gz", sample=SAMPLES),
		tempdir + "Probes/probes.sorted.bed.gz",
		tempdir + "Probes/probes.sorted.bed.gz.tbi",
		expand(tempdir + "ReadCounts/{sample}.num_reads.txt", sample=SAMPLES),
		tempdir + "TotalReads/total_read.txt",
		expand(tempdir + "RPM/{sample}.probe.rpm_rate.bed.gz", sample=SAMPLES),
		expand(tempdir + "RPM/{sample}.probe.rpm_rate.bed.gz.tbi", sample=SAMPLES),
		config['outputdir']

rule mosdepth:
	input:
		tempdir + "{sample}.bam"
	output:
		tempdir + "Mosdepth/{sample}.per-base.bed.gz",
		tempdir + "Mosdepth/{sample}.per-base.bed.gz.tbi"
	log:
		"logs/mosdepth.{sample}.log"
	params:
		tempdir = tempdir,
		outputdirname = tempdir + "Mosdepth/"
	shell:
		"""
		mkdir -p {params.outputdirname}
		mosdepth {params.outputdirname}{wildcards.sample} {input}
		tabix -p bed {params.outputdirname}{wildcards.sample}.per-base.bed.gz
		"""

rule count_reads:
	input:
		tempdir + "{sample}.bam"
	output:
		tempdir + "ReadCounts/{sample}.num_reads.txt"
	log:
		"logs/mosdepth.{sample}.log"
	params:
		tempdir = tempdir,
		outputdirname = tempdir + "ReadCounts/"
	shell:
		"""
		mkdir -p {params.outputdirname}
		samtools view -c -F 260 {input} > {output}
		"""

rule get_total_read_counts:
	input:
		expand(tempdir + "ReadCounts/{sample}.num_reads.txt", sample=SAMPLES)
	output:
		tempdir + "TotalReads/total_read.txt"
	params:
		tempdir = tempdir,
		outputdirname = tempdir + "TotalReads/"
	shell:
		"""
		mkdir -p {params.outputdirname}
		cat {input} > {params.outputdirname}total_read.txt
		"""

rule gzip_probes:
	input:
		probes = config['probes']
	output:
		tempdir + "Probes/probes.sorted.bed.gz",
		tempdir + "Probes/probes.sorted.bed.gz.tbi"
	params:
		tempdir = tempdir,
		outputdirname = tempdir + "Probes/"
	shell:
		"""
		mkdir -p {params.outputdirname}
		bedtools sort -i {input.probes} > {params.outputdirname}probes.sorted.bed
		bgzip {params.outputdirname}probes.sorted.bed
		tabix {params.outputdirname}probes.sorted.bed.gz -p bed
		"""

rule get_probe_reads_per_million:
	input:
		per_base_bed_gz=tempdir + "Mosdepth/{sample}.per-base.bed.gz",
		total_reads=tempdir + "TotalReads/total_read.txt",
		probes_gz=tempdir + "Probes/probes.sorted.bed.gz"
	output:
		tempdir + "RPM/{sample}.probe.rpm_rate.bed.gz",
		tempdir + "RPM/{sample}.probe.rpm_rate.bed.gz.tbi"
	params:
		tempdir = tempdir,
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

rule get_probe_cover_mean_std_for_reference_panel:
	#       Calculate various stats about RPM coverage by comparing one reference panel sampel to the rest of the referene panel
	input:
		ref_rpm=expand(tempdir + 'RPM/{sample}.probe.rpm_rate.bed.gz', sample=SAMPLES),
		ref_rpm_tbi=expand(tempdir + 'RPM/{sample}.probe.rpm_rate.bed.gz.tbi', sample=SAMPLES),
		sample=tempdir + 'RPM/{sample}.probe.rpm_rate.bed.gz',
		ref_sample_tbi=tempdir + 'RPM/{sample}.probe.rpm_rate.bed.gz.tbi'
	output:
		bedgz=tempdir + "ProbeCoverage/{sample}.probe.cover.mean.stdev.bed.gz",
		tbi=tempdir + "ProbeCoverage/{sample}.probe.cover.mean.stdev.bed.gz.tbi",
		bed=tempdir + "ProbeCoverage/{sample}.probe.cover.mean.stdev.bed"
	threads: 2
	params:
		tempdir = tempdir,
		outputdirname = tempdir + "ProbeCoverage/"
	shell:
		"""
		mkdir -p {params.outputdirname}
		get_regions_zscores.py -r {params.tempdir}RPM/  -s {input.sample} | bgzip -c > {output.bedgz}
		cp {output.bedgz} {params.outputdirname}{wildcards.sample}.probe.cover.mean.stdev.copy.bed.gz
		gunzip {params.outputdirname}{wildcards.sample}.probe.cover.mean.stdev.copy.bed.gz
		mv {params.outputdirname}{wildcards.sample}.probe.cover.mean.stdev.copy.bed {output.bed}
		tabix -p bed {output.bedgz}
		"""

rule get_adj_zscore_for_ref_panel:
	#       Calculate the z-score of the RPM data comparing each of the reference panel samples to the panel
	input:
		ref_sample_rpm=tempdir + "RPM/{sample}.probe.rpm_rate.bed.gz",
		ref_sample_coverage=tempdir + "ProbeCoverage/{sample}.probe.cover.mean.stdev.bed.gz",
		ref_sample_coverage_tbi=tempdir + "ProbeCoverage/{sample}.probe.cover.mean.stdev.bed.gz.tbi"
	output:
		bedgz=tempdir + "AdjZscore/{sample}.adj_z.bed.gz",
		tbi=tempdir + "AdjZscore/{sample}.adj_z.bed.gz.tbi"
	params:
		tempdir = tempdir,
		outputdirname = tempdir + "AdjZscore/"
	shell:
		"""
		mkdir -p {params.outputdirname}
		get_coverage_zscores.py \
			-r {input.ref_sample_rpm} \
			-s {input.ref_sample_coverage} \
			| bgzip -c > {output.bedgz}
		tabix -p bed {output.bedgz}
		"""

rule create_panel_db:
	input:
		rpm_gz=expand(tempdir + "RPM/{sample}.probe.rpm_rate.bed.gz", sample=SAMPLES),
		rpm_tbi=expand(tempdir + "RPM/{sample}.probe.rpm_rate.bed.gz.tbi", sample=SAMPLES),
		z_bedgz=expand(tempdir + "AdjZscore/{sample}.adj_z.bed.gz", sample=SAMPLES),
		z_tbi=expand(tempdir + "AdjZscore/{sample}.adj_z.bed.gz.tbi", sample=SAMPLES),
		all_calls=config['all_calls'],
		sample_list=SAMPLES_FILE
	output:
		directory(config['outputdir'])
	shell:
		"""
		mkdir -p {output}/RPM
		mkdir -p {output}/Calls
		mkdir -p {output}/AdjZscore
		cp {input.rpm_gz} {output}/RPM
		cp {input.rpm_tbi} {output}/RPM
		cp {input.z_bedgz} {output}/AdjZscore
		cp {input.z_tbi} {output}/AdjZscore	
		cp {input.all_calls} {output}/Calls
		cp {input.sample_list} {output}
		"""
