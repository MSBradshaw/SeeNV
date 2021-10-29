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

```
TODO
```

## Generate plots

```
cnviz \
-i cohort.samples \
-p probes.bed \
-c all_calls.txt \
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


