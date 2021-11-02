<p align="center"><img src="https://github.com/MSBradshaw/CNViz/blob/main/Logo.png?raw=true" width="80%"/></p>

CNViz is still being developed. It is can be download and used but is by no means exaustively tested. 

Find a bug? Create an [issue](https://github.com/MSBradshaw/CNViz/issues)! 

Have a feature idea? Create an [issue](https://github.com/MSBradshaw/CNViz/issues)!

CNViz provides comprehensive yet easy to digest visualizations for each call in a sample and depicts relevant statistics comparing a sample to a cohort of other samples. 
It is known that the accuracy and reliability of CNV calls increases when using multiple callers and parameter sets, for this reason CNViz also provides a way to visualize multiple callers or bin sizes simultaneously — a feature not known to exist in other tools. 
Combined with the tool [PlotCritic](https://github.com/jbelyeu/PlotCritic), we found that a clinician can accurately filter through roughly 200 calls in 20 minutes, or just 6 seconds per call on average. 
CNViz has been packaged as a Nextflow workflow that can be run as it’s own individual pipeline or included as part of another. 


<p align="center"><img src="https://github.com/MSBradshaw/CNViz/blob/main/CNVizPoster.png?raw=true" width="100%"/></p>


Poster Presentation from Genome Informatics 2021

# Installation

CNViz is currently only usable on Linux based systems. 

CNViz requires you have [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) installed.

Download this repo:

`git install git@github.com:MSBradshaw/CNViz.git`

Move into the repo:

`cd CNViz`

Run the install script:

`source install.sh`

The install script will create a conda envrionment called `cnviz` with all the necessary python dependancies for CNViz. It will also download a the following external non-python tools:

[Snakemake](https://snakemake.readthedocs.io/en/stable/index.html)

[gargs](https://github.com/brentp/gargs)

[mosdepth](https://github.com/brentp/mosdepth)

[samtools](http://www.htslib.org/)

[bedtools](https://bedtools.readthedocs.io/en/latest/)

All these dependacies will be placed in the bin direcoty of the cnviz conda environment.

As long as the conda environemnt is activated the `cnviz` command can now be used anywhere.

# Usage

CNViz requires it's conda environment to work, start the conda environment:

`conda activate cnviz`

## Build a Reference DB

In order to generate plots, a reference panel database is required. You can either create your own or use the one included with this repository.

```
cnviz \
-b \
-i panel.samples \
-s ../WES_TargetCoverage_v2.bed \
-c ../all_calls.txt \
-o DELDELDEL \
-t 32
```

## Generate plots

```
cnviz \
-p \
-i Example/cohort.samples \
-s Example/probes.bed \
-c Example/all_calls.txt \
-o OutputDir \
-r RefPanel \
-g Homo_sapiens.GRCh37.82.genes.bed.gz \
-a gnomad_v2.1_sv.sites.bed.gz \
-t 50
```

## Parameter explaination
```
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
```

## `-i INPUT_SAMPLES`

Same format for samples and reference panel. See `Example/panel.samples` or `Example/cohort.samples`

TSV file with the following columns in this order:

Patient Id: Unique identifier for each patient

Sample Id: Unique identifier for each sample (can be the same as the sample's Patient Id)

Sex: Chromosomal sex of patient (e.g. XY or XX). If unknown just input ZZ.

BAM: path to the sample's .bam file

BAI: path to the sample's .bai file

VCF: path to the sample's .vcf file

TBI: path to the sample's .vcf.tbi file

## `-s SITES`

[BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) formatted file with sites of interest within the genome. We use the regions our probes target for whole exome sequencing. See `Example/probes.bed`

## `-c ALL_CALLS `

File that contains a list of files containing calls. See `Example/all_calls.txt`.

Each file contains the path to a set of calls in a bed file with the following information seporated by tabs:

Chromosome: the Chromosome number/name (if listing chromsome 1 input 1 not chr1)

Start: base number at which the calls starts

End: base number at which the calls ends

Call type: type of call made (Duplication, Deletion ect)

Sample Id: Sample Id the call pertains to, should match Sample Id's found in `-i INPUT_SAMPLES` 

## `-o OUTPUT`

When using the `-b` or `--buildref` flag `-o` will is the path to where you want the reference panel database (a dirrectory) saved.

When using the `-p` or `--plotsamples` flag `-o` is the path for where the plots will be saved.

## `-r REF_DB`

Path to a reference panel database (the output created when using the `-b` flag)

## `-g GENES_FILE`

Path to a gzipped or bgzipped bed formated file file with the columns:

Chromosome: the Chromosome number/name (if listing chromsome 1 input 1 not chr1)

Start: base number at which the calls starts

End: base number at which the calls ends

Gene symbol: the gene symbol

You can either download data and format it from [biomart](http://uswest.ensembl.org/biomart/martview/) or use our HG19 based verion found in Example/hg19.genes.bed.gz

## `-a GNOMAD with Allele Frequency data`

path to the gnomAD SV sites file in .beg.gz format with an accompanying .beg.gz.tbi file in the same dirrectory. Download [.beg.gz file](https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.bed.gz) and [.beg.gz.tbi](https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.bed.gz.tbi)

## `-t THREADS`

The number of _cores_ (not threads) to be used. Default is 1, but you really should use many more. For reference, using 32 cores, the provided reference panel database, and processing/plotting 6 samples with 300 calls takes us ~2 hours.
