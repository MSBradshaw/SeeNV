import pandas as pd

samples_df = pd.read_csv('panel.samples',sep='\t',header=None)
samples_df.columns = ['p_id','s_id','sex','bam','bai','vcf','tbi']

SAMPLES = samples_df['s_id']
BAM = samples_df["bam"]
print(SAMPLES)
rule all:
	input:
		expand("{sample}.bam", sample=SAMPLES)

rule work_with_bams:
	input:
		expand("{bam}", zip, bam=BAM, sample)
	params:
		sid=expand("{sample}.bam {bam}", zip, sample=SAMPLES, bam=BAM)
	output:
		"{sample}.bam"
	shell:
		"ln -s {sid} {output}"
