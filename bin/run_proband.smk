import pandas as pd
import os

# Load input files and options
samples_df = pd.read_csv('workproband/proband.samples',sep='\t',header=None)
samples_df.columns = ['p_id','s_id','sex','bam','bai']

SAMPLES = samples_df['s_id']
BAM = samples_df["bam"]

for line in open('workproband/outputdir.txt','r'):
	config['outputdir'] = line.strip()
	break
for line in open('workproband/reference_panel.txt','r'):
	config['reference_panel'] = line.strip()
	if config['reference_panel'][-1] != '/':
		config['reference_panel'] += '/'
	break

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


rule all:
	input:
		expand("workproband/Mosdepth/{sample}.per-base.bed.gz", sample=SAMPLES),
		"workproband/Probes/probes.sorted.bed.gz",
                "workproband/Probes/probes.sorted.bed.gz.tbi",
		expand("workproband/ReadCounts/{sample}.num_reads.txt", sample=SAMPLES),
		"workproband/TotalReads/total_read.txt",
		expand("workproband/RPM/{sample}.probe.rpm_rate.bed.gz", sample=SAMPLES),
                expand("workproband/RPM/{sample}.probe.rpm_rate.bed.gz.tbi", sample=SAMPLES),
		"workproband/Params/call_plotting_params.txt",
		#"workproband/Exons/labeled_exons.bed.gz",
		#"workproband/Exons/labeled_exons.bed.gz.tbi",
		expand("workproband/ProbeCoverage/{sample}_probe.cover.mean.stdev.bed.gz", sample=SAMPLES),
		expand("workproband/ProbeCoverage/{sample}_probe.cover.mean.stdev.bed.gz.tbi", sample=SAMPLES),
		expand('workproband/ProbeCoverageRefPanel/{ref_sample}.probe.rpm_rate.bed.gz', ref_sample=REF_SAMPLES),
                expand('workproband/ProbeCoverageRefPanel/{ref_sample}.probe.rpm_rate.bed.gz.tbi', ref_sample=REF_SAMPLES),
		expand("workproband/AdjZscore/{sample}.adj_z.bed.gz", sample=SAMPLES),
		expand("workproband/AdjZscore/{sample}.adj_z.bed.gz.tbi", sample=SAMPLES),
		"workproband/MergedAdjZscore/adj_scores.bed.gz",
		"workproband/MergedAdjZscore/adj_scores.bed.gz.tbi",
		"workproband/SplitParams/.split_param_file.ready",
		expand("workproband/PlotsComplete/{num}.done",num=NUM_CALLS),
		expand("workproband/AdjZscoreRefPanel/{ref_sample}.adj_z.bed.gz", ref_sample=REF_SAMPLES),
		expand("workproband/AdjZscoreRefPanel/{ref_sample}.adj_z.bed.gz.tbi", ref_sample=REF_SAMPLES),
		expand("workproband/AlleleCount/{sample}.allele_count.tsv", sample=SAMPLES)


rule mosdepth:
	input:
		"workproband/{sample}.bam"
	output:
		"workproband/Mosdepth/{sample}.per-base.bed.gz",
		"workproband/Mosdepth/{sample}.per-base.bed.gz.tbi"
	log:
		"workproband/logs/mosdepth.{sample}.log"
	shell:
		"""
		mkdir -p workproband/Mosdepth
		mosdepth workproband/Mosdepth/{wildcards.sample} {input}
		tabix -f -p bed workproband/Mosdepth/{wildcards.sample}.per-base.bed.gz
		"""

rule count_reads:
#	"""
#	Count the number of reads in each sample's bam
#	"""
        input:
                "workproband/{sample}.bam"
        output:
                "workproband/ReadCounts/{sample}.num_reads.txt"
        log:
                "workproband/logs/count_reads.{sample}.log"
        shell:
            """
			mkdir -p workproband/ReadCounts
			samtools view -c -F 260 {input} > {output}
            """

rule get_total_read_counts:
#	"""
#	Add all the individual read counts to one file	
#	"""
	input:
		expand("workproband/ReadCounts/{sample}.num_reads.txt", sample=SAMPLES)
	output:
		"workproband/TotalReads/total_read.txt"
	shell:
		"""
		mkdir -p workproband/TotalReads
		cat {input} > workproband/TotalReads/total_read.txt
		"""

rule gzip_probes:
#	"""
#	Take the probes files, bgzip and tabix it
#	"""
	input:
		probes = config['probes']
	output:
		"workproband/Probes/probes.sorted.bed.gz",
		"workproband/Probes/probes.sorted.bed.gz.tbi"
	shell:
		"""
		mkdir -p workproband/Probes/
		bedtools sort -i {input.probes} > workproband/Probes/probes.sorted.bed
		bgzip workproband/Probes/probes.sorted.bed
		tabix workproband/Probes/probes.sorted.bed.gz -p bed
		"""

rule get_probe_reads_per_million:
#	"""
#	For each sample calculate the number of reads per million (RPM) at each probe. Return the files bgzip-ed and tabix-ed
#	"""
	input:
		per_base_bed_gz="workproband/Mosdepth/{sample}.per-base.bed.gz",
		total_reads="workproband/TotalReads/total_read.txt",
		probes_gz="workproband/Probes/probes.sorted.bed.gz"
	output:
		"workproband/RPM/{sample}.probe.rpm_rate.bed.gz",
		"workproband/RPM/{sample}.probe.rpm_rate.bed.gz.tbi"
	shell:
		"""
		mkdir -p workproband/RPM/
		get_rpm_rates.py \
		-m {input.per_base_bed_gz} \
		--regions_file {input.probes_gz} \
		--num_reads {input.total_reads} \
		| bgzip -c > workproband/RPM/{wildcards.sample}.probe.rpm_rate.bed.gz
		tabix -p bed workproband/RPM/{wildcards.sample}.probe.rpm_rate.bed.gz
		"""

rule create_plotting_params:
#	"""
#	Create a file will add the CNVs and their preconficgures parameters needed for plotting
#	"""
	input:
		all_calls=config['all_calls']

	output:
		"workproband/Params/call_plotting_params.txt"
	params:
		samples=expand("{sample}",sample=SAMPLES)
	log:
		"workproband/logs/plotting_params.log"
	shell:
		"""
		mkdir -p workproband/Params/
		cat {input.all_calls} | awk '{{ print $1":"$2"-"$3"\t"$4"\t"$5"\t.\t"$1"."$2"-"$3"\t"$4}}' >> workproband/Params/tmp_call_plotting_params.txt
		
		for i in {{1..23}}; do
			grep "^$i" workproband/Params/tmp_call_plotting_params.txt >> workproband/Params/call_plotting_params.txt.tmp || echo "Not found" 2>> {log}
		done

		touch workproband/Params/call_plotting_params.txt
		# remove samples from the calls that are not in the sample list
		array=( {params.samples} )
		for j in "${{array[@]}}"
		do
			cat workproband/Params/call_plotting_params.txt.tmp | {{ grep $j || :; }} >> workproband/Params/call_plotting_params.txt
		done
		"""

checkpoint split_param_file:
#	"""
#	Take the one large param files and split it up into files with a single line in it
#	"""
	input:
		param_file="workproband/Params/call_plotting_params.txt",
	output:
		"workproband/SplitParams/.split_param_file.ready"
	params:
		samples=expand("{sample}",sample=SAMPLES)
	log:
		"workproband/logs/split_param_file.log"
	shell:
		"""
		mkdir -p workproband/SplitParams
		i=0
		while read p; do
			echo $p > workproband/SplitParams/${{i}}.txt
			i=$((i+1))
		done <{input.param_file}
		touch workproband/SplitParams/.split_param_file.ready
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

		names = glob_wildcards('workproband/SplitParams/{rs}.txt')[0]

		return self.pattern

rule label_exons:
	"""
	Take in the genes files and add create named labels for the exons in all genes.
	Return it bgzip-ed and tabix-ed
	"""
	input:
		"workproband/genes.bed.gz"
	output:
		"workproband/Exons/labeled_exons.bed.gz",
		"workproband/Exons/labeled_exons.bed.gz.tbi"
	shell:
		"""
		mkdir -p workproband/Exons/
		zcat {input} | label_exon_number.py > workproband/Exons/labeled_exons.bed
		bgzip workproband/Exons/labeled_exons.bed
		tabix -p bed workproband/Exons/labeled_exons.bed.gz
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
		'workproband/ProbeCoverageRefPanel/{ref_sample}.probe.rpm_rate.bed.gz',
		'workproband/ProbeCoverageRefPanel/{ref_sample}.probe.rpm_rate.bed.gz.tbi',
		'workproband/AdjZscoreRefPanel/{ref_sample}.adj_z.bed.gz',
		'workproband/AdjZscoreRefPanel/{ref_sample}.adj_z.bed.gz.tbi'
	shell:
		"""
		mkdir -p workproband/ProbeCoverageRefPanel/
		mkdir -p workproband/AdjZscoreRefPanel/
		cp {input.bedgz} workproband/ProbeCoverageRefPanel/
		cp {input.tbi} workproband/ProbeCoverageRefPanel/
		cp {input.adj_bed} workproband/AdjZscoreRefPanel/
		cp {input.adj_tbi} workproband/AdjZscoreRefPanel/
		"""


rule get_probe_cover_mean_std:
	"""
	Calculate various stats about RPM coverage by comparing one proband to the reference panel
	"""
	input:
		ref_rpm=expand('workproband/ProbeCoverageRefPanel/{ref_sample}.probe.rpm_rate.bed.gz', ref_sample=REF_SAMPLES),
		ref_rpm_tbi=expand('workproband/ProbeCoverageRefPanel/{ref_sample}.probe.rpm_rate.bed.gz.tbi', ref_sample=REF_SAMPLES),
		proband_rpm="workproband/RPM/{sample}.probe.rpm_rate.bed.gz"	
	output:
		bedgz="workproband/ProbeCoverage/{sample}_probe.cover.mean.stdev.bed.gz",
		tbi="workproband/ProbeCoverage/{sample}_probe.cover.mean.stdev.bed.gz.tbi",
		bed="workproband/ProbeCoverage/{sample}_probe.cover.mean.stdev.bed"
	threads: 2
	shell:
		"""
		mkdir -p workproband/ProbeCoverage/
		get_regions_zscores.py -r workproband/ProbeCoverageRefPanel/ -s {input.proband_rpm} | bgzip -c > {output.bedgz}
		cp {output.bedgz} workproband/ProbeCoverage/{wildcards.sample}_probe.cover.mean.stdev.copy.bed.gz
		gunzip workproband/ProbeCoverage/{wildcards.sample}_probe.cover.mean.stdev.copy.bed.gz
		mv workproband/ProbeCoverage/{wildcards.sample}_probe.cover.mean.stdev.copy.bed {output.bed}
		tabix -p bed {output.bedgz}
		"""


rule get_adj_zscore:
#	"""
#	Calculate the z-score of the RPM data compared to the std and mean for each proband 
#	"""
	input:
		proband_rpm="workproband/RPM/{sample}.probe.rpm_rate.bed.gz",
		proband_coverage="workproband/ProbeCoverage/{sample}_probe.cover.mean.stdev.bed.gz",
		proband_coverage_tbi="workproband/ProbeCoverage/{sample}_probe.cover.mean.stdev.bed.gz.tbi"
	output:
		bedgz="workproband/AdjZscore/{sample}.adj_z.bed.gz",
		tbi="workproband/AdjZscore/{sample}.adj_z.bed.gz.tbi"
	shell:
		"""
		mkdir -p workproband/AdjZscore/
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
		bedgz=expand("workproband/AdjZscore/{sample}.adj_z.bed.gz", sample=SAMPLES),
		tbi=expand("workproband/AdjZscore/{sample}.adj_z.bed.gz.tbi", sample=SAMPLES),
		ref_bed=expand("workproband/AdjZscoreRefPanel/{ref_sample}.adj_z.bed.gz", ref_sample=REF_SAMPLES),
		ref_tbi=expand("workproband/AdjZscoreRefPanel/{ref_sample}.adj_z.bed.gz.tbi", ref_sample=REF_SAMPLES)
	output:
		final_bed_gz="workproband/MergedAdjZscore/adj_scores.bed.gz",
		final_bed_gz_tbi="workproband/MergedAdjZscore/adj_scores.bed.gz.tbi"
	shell:
		"""
		mkdir -p workproband/AllAdjZScore
		cp workproband/AdjZscore/* workproband/AllAdjZScore
		cp workproband/AdjZscoreRefPanel/* workproband/AllAdjZScore
		merge_samples_adj_scores.py -r workproband/AllAdjZScore | bgzip -c > {output.final_bed_gz}
		tabix -p bed {output.final_bed_gz}
		"""

rule prep_all_calls:
	input:
		all_calls=config['all_calls']
	output:
		"workproband/FindMaxTMP/.ready"
	params:	
		samples=expand("{sample}",sample=SAMPLES)
	log:
		"workproband/logs/prep_all_calls.log"
	shell:
		"""
		mkdir -p workproband/FindMaxTMP
		rm -f workproband/FindMaxTMP/*
		
		name=$(basename {input.all_calls})
		array=( {params.samples} )
		for i in "${{array[@]}}"
		do
			cat {input.all_calls} | {{ grep $i || :; }} >> workproband/FindMaxTMP/$name.filtered 2>> {log}
		done

		bedtools sort -i workproband/FindMaxTMP/$name.filtered > workproband/FindMaxTMP/$name.sorted 2>> {log}
		bgzip -f workproband/FindMaxTMP/$name.sorted 2>> {log}
		tabix -f -p bed workproband/FindMaxTMP/$name.sorted.gz 2>> {log}
		echo "workproband/FindMaxTMP/${{name}}.sorted.gz" >> workproband/FindMaxTMP/multiple_savvy_calls.txt 2>> {log}
		
		touch workproband/FindMaxTMP/.ready
		"""

rule plotter:
	input:
		probes="workproband/Probes/probes.sorted.bed.gz",
		probes_tbi="workproband/Probes/probes.sorted.bed.gz.tbi",
		all_calls="workproband/FindMaxTMP/.ready",
		param_file=Checkpoint_MakePattern_prime('workproband/SplitParams/{num}.txt'),
		marker_file="workproband/SplitParams/.split_param_file.ready",
		coverage=expand("workproband/ProbeCoverage/{sample}_probe.cover.mean.stdev.bed.gz", sample=SAMPLES),
		coverage_tbi=expand("workproband/ProbeCoverage/{sample}_probe.cover.mean.stdev.bed.gz.tbi", sample=SAMPLES),
		adj_scores="workproband/MergedAdjZscore/adj_scores.bed.gz",
		adj_scores_tbi="workproband/MergedAdjZscore/adj_scores.bed.gz.tbi",
		gnomad_sv=config['gnomad_sv_bed_gz'],
		vardb=config['vardb_bed_gz'],
		repeat_masker=config['repeat_masker_gz']

	output:
		"workproband/PlotsComplete/{num}.done"
	log:
		"workproband/logs/plotting_{num}.log"
	params:
		num_calls=len(NUM_CALLS),
		site_quality=config['site_quality']
	shell:
		"""
		mkdir -p workproband/PlotsComplete/
		inputs=""
		empty=""
		while read p; do
			if [ "$inputs" == "$empty" ]; then
				inputs="$p"
			else
				inputs="$inputs,$p"
			fi
		done < workproband/FindMaxTMP/multiple_savvy_calls.txt
		echo {params.num_calls} >> {log}
		echo {wildcards.num} >> {log}
		mkdir -p workproband/Plots/ 2>> {log}
		mkdir -p workproband/PlotsComplete/ 2>> {log}
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
			--output workproband/Plots/"${{sample}}.${{region}}.${{svtype}}".png \
			--coverage_scores {input.adj_scores} \
			--sites {input.probes} \
			--window 50000 \
			--region ${{region}} \
			--title "${{sample}} ${{region}} ${{svtype}}" \
			--depth workproband/ProbeCoverage/${{sample}}_probe.cover.mean.stdev.bed \
			--gnomad {input.gnomad_sv} \
			--vardb {input.vardb} \
			--repeatmasker {input.repeat_masker} \
			--site_quality {params.site_quality}
		touch {output}
		"""

rule collect_allele_counts:
	input:
		bam="workproband/{sample}.bam",
		probes=config['probes']
	output:
		"workproband/AlleleCount/{sample}.allele_count.tsv"
	resources:
		mem_mb=20000
	shell: "mkdir -p workproband/AlleleCount; touch {output}"
