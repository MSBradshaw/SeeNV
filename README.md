# CNViz
<img src="https://github.com/MSBradshaw/CNViz/blob/main/Logo.png?raw=true" width="600"/>

CNViz provides comprehensive yet easy to digest visualizations for each call in a sample and depicts relevant statistics comparing a sample to a cohort of other samples. 
It is known that the accuracy and reliability of CNV calls increases when using multiple callers and parameter sets, for this reason CNViz also provides a way to visualize multiple callers or bin sizes simultaneously — a feature not known to exist in other tools. 
Combined with the tool [PlotCritic](https://github.com/jbelyeu/PlotCritic), we found that a clinician can accurately filter through roughly 200 calls in 20 minutes, or just 6 seconds per call on average. 
CNViz has been packaged as a Nextflow workflow that can be run as it’s own individual pipeline or included as part of another. 

<img src="https://github.com/MSBradshaw/CNViz/blob/main/CNVizPoster.png?raw=true" width="600"/>

Poster Presentation from Genome Informatics 2021

## Installation

CNViz depends on just two things

1. [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
2. [Singularity](https://sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps)

All other dependancies and packages are contained in a single singularity container.

Clone this repo

`git clone git@github.com:MSBradshaw/CNViz.git`

Build singularity container (this will require sudo privlidges)

`cd CNViz/Singularity`

`bash build.sh`

or if you do not have sudo privlidges, download a prebuilt version of the container from [here](https://drive.google.com/file/d/1qaX7MfytQttYyVetnVdtoHtfx7woTpf4/view?usp=sharing):

## Configuration

Nextflow requires some configuration. This part of the documentation is forth coming. If you happen to be experienced with nextflow already, AWESOME! You probbaly already know what to do!

## Usuage

In order to use CNViz a reference panel DB is required. You can either use your own, or use this example one provided here (sorry the example I am legally allowed to publicly share is forth coming):

### Build your own reference panel DB

TODO

### Run CNViz pipeline for probands

TODO

TODO: pretty example output pictures



<img src="https://github.com/MSBradshaw/CNViz/blob/main/CNViz.png?raw=true" width="600"/>
