<p align="center"><img src="https://github.com/MSBradshaw/SeeNV/blob/main/Imgs/SeeNVLogoBlue.png?raw=true" width="80%"/></p>

<p align="center"><img src="https://github.com/MSBradshaw/SeeNV/blob/main/Imgs/dup.png?raw=true" width="80%"/></p>

SeeNV is a tool for visualizing and assessing the technical evidence behind a copy number variation (CNV) call identified in whole exome sequencing data.
For each CNV in the input, SeeNV generates an infographic with statistics about the normalized coverage, controlling for variation within individual samples, a reference database of samples, and the sequencing batch. Additionally, it includes information about the abundance of the call in the wider population and information about the complexity and variability of coverage in that region of the genome.
Using SeeNV you can rapidly and reliably assess the validitity of CNV calls, we found onaverage a user needs only 7.5 seconds to determine if a call is real based on the infographs and achieves 0.93 precision and 0.72 recall.

Find a bug? Create an [issue](https://github.com/MSBradshaw/SeeNV/issues)! 

Have a feature idea? Create an [issue](https://github.com/MSBradshaw/SeeNV/issues)!

# Installation

SeeNV is currently only usable on Linux-based systems. 

SeeNV requires you to have [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) installed.

Download this repo, move into it, and run the installation script!

```
git clone https://github.com/MSBradshaw/SeeNV.git
cd SeeNV
source install.sh
```

The install script will create a conda environment called `seenv` with all the necessary python dependencies for SeeNV. It will also download the following external non-python tools:

[Snakemake](https://snakemake.readthedocs.io/en/stable/index.html)

[gargs](https://github.com/brentp/gargs)

[mosdepth](https://github.com/brentp/mosdepth)

[samtools](http://www.htslib.org/)

[bedtools](https://bedtools.readthedocs.io/en/latest/)

All these dependencies will be placed in the bin directory of the seenv conda environment.

As long as the conda environment is activated the `seenv` command can now be used anywhere.

# Usage

SeeNV requires its conda environment to work, start the conda environment:

`conda activate seenv`

## Build a Reference DB

In order to generate plots, a reference panel database is required. You can either create your own or use the one included in this repository. 

Note: If running the example besure to download the 1000 Genomes Project samples using `ExampleData/download_1KG_samples.sh` ([aws cli](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html) required) and change the paths in `ExampleData/ref_panel_sample_list.tsv` to point to the proper locations.

```
seenv \
-b \
-i ExampleData/ref_panel_sample_list.tsv \
-s ExampleData/SureSelect_All_Exon_V2.bed \
-c ExampleData/17_bi_300_samples_calls.bed \
-o ReferenceDB \
-t 32
```

## Generate plots

```
seenv \
-p \
-i ExampleData/proband_sample_list.tsv \
-s ExampleData/SureSelect_All_Exon_V2.bed \
-c ExampleData/17_bi_300_samples_calls.bed \
-a ExampleData/gnomad_v2.1_sv.sites.bed.gz \
-t 64 \
-v ExampleData/vardb.sorted.bed.gz \
-m ExampleData/genomicRepeats.sorted.bed.gz \
-r ReferenceDB/ \
-o TestPlots
```

## Parameter explanation
```
One of the following options is required

--buildref, -b      flag to use function to build a reference panel database
--plotsamples, -p   flag to plot samples

Parameters to accompany --plotsamples, -p:
    -i INPUT_SAMPLES   (required) samples list
    -s SITES           (required) genomic sites bed file
    -c CALLS           (required) bed file containing all calls with columns: chrom, start, end, cnv_type, sample_name
    -o OUTPUT          (required) output directory, where to save the plots
    -r REF_DB          (required) path to reference db created by the --buildref function of SeeNV
    -a GNOMAD          (required) the gnomad sv sites file with allele frequency information
    -t THREADS         (optional) number of threads to use, default 1 (you really want to use more than 1)
    -v varDB           (required) path to a GZipped bed file for the varDB common variants with an accompanying tabix indexed
    -m RepeatMasker    (required) path to a GZipped bed file for the RepeatMasker elements with an accompanying tabix indexed
    -q Site Quality    (optional) index of the column in the sites file that contains the site quality statement, default None

Parameters to accompany --buildref, -b
    -i INPUT_SAMPLES  (required) samples list
    -s SITES          (required) genomic sites bed file
    -c CALLS          (required) bed file containing all calls with columns: chrom, start, end, cnv_type, sample_name
    -o OUTPUT         (required) output directory, reference panel database
    -t THREADS        (optional) number of threads to use, default 1 (you really want to use more than 1)
```

## `-i INPUT_SAMPLES`

Same format for samples and reference panel. See `Example/ref_sample_list.tsv` or `Example/proband_sample_list.tsv`

TSV file with the following columns in this order:

Patient Id: Unique identifier for each patient

Sample Id: Unique identifier for each sample (can be the same as the sample's Patient Id)

Sex: Chromosomal sex of patient (e.g. XY or XX). If unknown just input ZZ.

BAM: absolute path to the sample's .bam file

BAI: absolute path to the sample's .bai file

## `-s SITES`

[BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) formatted file with sites of interest within the genome. We use the regions our probes target for whole exome sequencing. See `Example/probes.bed`

## -q Site Quality    

(optional) the index of the column in the sites file that contains the site quality statement, default None. You can add an additional column to your sites bed file, in this column you can label probes as "Good" (case insensitive) if it is a probe/site that passes your own quality control metrics. This information will be used to populate the probe quality bar under the main plots. Non "Good" strings are treated equivalently as not-good.

## `-c CALLS `

See example `ExampleData/17_bi_300_samples_calls.bed`

Bed file with the following information separated by tabs:

Chromosome: the Chromosome number/name (if listing chromosome 1 input 1, not chr1)

Start: base number at which the call starts

End: base number at which the call ends

Call type: type of call made (Duplication, Deletion ect)

Sample Id: Sample Id the call pertains to, should match Sample Ids found in `-i INPUT_SAMPLES` 

## `-o OUTPUT`

When using the `-b` or `--buildref` flag `-o` will is the path to where you want the reference panel database (a dirrectory) saved.

When using the `-p` or `--plotsamples` flag `-o` is the path for where the plots will be saved.

## `-r REF_DB`

Path to a reference panel database (the output created when using the `-b` flag)

## `-a GNOMAD with Allele Frequency data`

path to the gnomAD SV sites file in .beg.gz format with an accompanying .beg.gz.tbi file in the same directory. This can be found in `ExampleData/gnomad_v2.1_sv.sites.bed.gz` or can the `bed.gz` and and `.tbi` be downloaded from [here](https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.bed.gz) and [here](https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.bed.gz.tbi)

## `-m RepeatMasker `

path to a GZipped bed file for the varDB common variants with an accompanying tabix indexed. This can be found in `ExampleData/genomicRepeats.sorted.bed.gz`.

## `-v varDB`

path to a GZipped bed file for the varDB common variants with an accompanying tabix indexed. This can be found in `ExampleData/vardb.sorted.bed.gz`

## `-t THREADS`

The number of _cores_ (not threads) to be used. The default is 1, but you really should use many more. For reference, using 32 cores, the provided reference panel database, and processing/plotting 6 samples with 300 calls takes us ~2 hours.

# Example output figures
## Deletion
<p align="center"><img src="https://github.com/MSBradshaw/SeeNV/blob/main/Imgs/del.png?raw=true" width="80%"/></p>

## Duplication
<p align="center"><img src="https://github.com/MSBradshaw/SeeNV/blob/main/Imgs/dup.png?raw=true" width="80%"/></p>
