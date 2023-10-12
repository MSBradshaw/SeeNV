import pandas as pd

SAMPLES_FILE = "workpanel/panel.samples"

# Load input files and options
samples_df = pd.read_csv(SAMPLES_FILE,sep='\t',header=None)
samples_df.columns = ['p_id','s_id','sex','bam','bai']

PROBESFILE = 'workpanel/Probes/probes.original.bed'

SAMPLES = samples_df['s_id']
BAM = samples_df["bam"]

for line in open('workpanel/outputdir.txt','r'):
	config['outputdir'] = line.strip()
	break

print(config)

rule all:
	input:
		expand("workpanel/Mosdepth/{sample}.per-base.bed.gz", sample=SAMPLES),
		"workpanel/Probes/probes.sorted.bed.gz",
                "workpanel/Probes/probes.sorted.bed.gz.tbi",
		expand("workpanel/ReadCounts/{sample}.num_reads.txt", sample=SAMPLES),
		"workpanel/TotalReads/total_read.txt",
		expand("workpanel/RPM/{sample}.probe.rpm_rate.bed.gz", sample=SAMPLES),
                expand("workpanel/RPM/{sample}.probe.rpm_rate.bed.gz.tbi", sample=SAMPLES),
		config['outputdir']

rule mosdepth:
	input:
		"workpanel/{sample}.bam"
	output:
		"workpanel/Mosdepth/{sample}.per-base.bed.gz",
		"workpanel/Mosdepth/{sample}.per-base.bed.gz.tbi"
	log:
		"logs/mosdepth.{sample}.log"
	shell:
		"""
		mkdir -p workpanel/Mosdepth
		mosdepth workpanel/Mosdepth/{wildcards.sample} {input}
		tabix -p bed workpanel/Mosdepth/{wildcards.sample}.per-base.bed.gz
		"""

rule count_reads:
        input:
                "workpanel/{sample}.bam"
        output:
                "workpanel/ReadCounts/{sample}.num_reads.txt"
        log:
                "logs/mosdepth.{sample}.log"
        shell:
			"""
			mkdir -p workpanel/ReadCounts
			samtools view -c -F 260 {input} > {output}
			"""

rule get_total_read_counts:
	input:
		expand("workpanel/ReadCounts/{sample}.num_reads.txt", sample=SAMPLES)
	output:
		"workpanel/TotalReads/total_read.txt"
	shell:
		"""
		mkdir -p workpanel/TotalReads
		cat {input} > workpanel/TotalReads/total_read.txt
		"""

rule gzip_probes:
	input:
		probes = config['probes']
	output:
		"workpanel/Probes/probes.sorted.bed.gz",
		"workpanel/Probes/probes.sorted.bed.gz.tbi"
	shell:
		"""
		mkdir -p workpanel/Probes/
		bedtools sort -i {input.probes} > workpanel/Probes/probes.sorted.bed
		bgzip workpanel/Probes/probes.sorted.bed
		tabix workpanel/Probes/probes.sorted.bed.gz -p bed
		"""

rule get_probe_reads_per_million:
	input:
		per_base_bed_gz="workpanel/Mosdepth/{sample}.per-base.bed.gz",
		total_reads="workpanel/TotalReads/total_read.txt",
		probes_gz="workpanel/Probes/probes.sorted.bed.gz"
	output:
		"workpanel/RPM/{sample}.probe.rpm_rate.bed.gz",
		"workpanel/RPM/{sample}.probe.rpm_rate.bed.gz.tbi"
	shell:
		"""
		mkdir -p workpanel/RPM
		get_rpm_rates.py \
		-m {input.per_base_bed_gz} \
		--regions_file {input.probes_gz} \
		--num_reads {input.total_reads} \
		| bgzip -c > workpanel/RPM/{wildcards.sample}.probe.rpm_rate.bed.gz
		tabix -p bed workpanel/RPM/{wildcards.sample}.probe.rpm_rate.bed.gz
		"""

rule get_probe_cover_mean_std_for_reference_panel:
	#       Calculate various stats about RPM coverage by comparing one reference panel sampel to the rest of the referene panel
	input:
		ref_rpm=expand('workpanel/RPM/{sample}.probe.rpm_rate.bed.gz', sample=SAMPLES),
		ref_rpm_tbi=expand('workpanel/RPM/{sample}.probe.rpm_rate.bed.gz.tbi', sample=SAMPLES),
		sample='workpanel/RPM/{sample}.probe.rpm_rate.bed.gz',
		ref_sample_tbi='workpanel/RPM/{sample}.probe.rpm_rate.bed.gz.tbi'
	output:
		bedgz="workpanel/ProbeCoverage/{sample}.probe.cover.mean.stdev.bed.gz",
		tbi="workpanel/ProbeCoverage/{sample}.probe.cover.mean.stdev.bed.gz.tbi",
		bed="workpanel/ProbeCoverage/{sample}.probe.cover.mean.stdev.bed"
	threads: 2
	shell:
		"""
		mkdir -p workpanel/ProbeCoverage
		get_regions_zscores.py -r workpanel/ProbeCoverage -s {input.sample} | bgzip -c > {output.bedgz}
		cp {output.bedgz} workpanel/ProbeCoverage/{wildcards.sample}.probe.cover.mean.stdev.copy.bed.gz
		gunzip workpanel/ProbeCoverage/{wildcards.sample}.probe.cover.mean.stdev.copy.bed.gz
		mv workpanel/ProbeCoverage/{wildcards.sample}.probe.cover.mean.stdev.copy.bed {output.bed}
		tabix -p bed {output.bedgz}
		"""

rule get_adj_zscore_for_ref_panel:
	#       Calculate the z-score of the RPM data comparing each of the reference panel samples to the panel
	input:
		ref_sample_rpm="workpanel/RPM/{sample}.probe.rpm_rate.bed.gz",
		ref_sample_coverage="workpanel/ProbeCoverage/{sample}.probe.cover.mean.stdev.bed.gz",
		ref_sample_coverage_tbi="workpanel/ProbeCoverage/{sample}.probe.cover.mean.stdev.bed.gz.tbi"
	output:
		bedgz="workpanel/AdjZscore/{sample}.adj_z.bed.gz",
		tbi="workpanel/AdjZscore/{sample}.adj_z.bed.gz.tbi"
	shell:
		"""
		mkdir -p workpanel/ProbeCoverage
		get_coverage_zscores.py \
			-r {input.ref_sample_rpm} \
			-s {input.ref_sample_coverage} \
			| bgzip -c > {output.bedgz}
		tabix -p bed {output.bedgz}
		"""

rule create_panel_db:
	input:
		rpm_gz=expand("workpanel/RPM/{sample}.probe.rpm_rate.bed.gz", sample=SAMPLES),
		rpm_tbi=expand("workpanel/RPM/{sample}.probe.rpm_rate.bed.gz.tbi", sample=SAMPLES),
		z_bedgz=expand("workpanel/AdjZscore/{sample}.adj_z.bed.gz", sample=SAMPLES),
		z_tbi=expand("workpanel/AdjZscore/{sample}.adj_z.bed.gz.tbi", sample=SAMPLES),
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
