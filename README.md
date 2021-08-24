# CNViz

CNViz provides comprehensive yet easy to digest visualizations for each call in a sample and depicts relevant statistics comparing a sample to a cohort of other samples. 
It is known that the accuracy and reliability of CNV calls increases when using multiple callers and parameter sets, for this reason CNViz also provides a way to visualize multiple callers or bin sizes simultaneously — a feature not known to exist in other tools. 
Combined with the tool [PlotCritic](https://github.com/jbelyeu/PlotCritic), we found that a clinician can accurately filter through roughly 200 calls in 20 minutes, or just 6 seconds per call on average. 
CNViz has been packaged as a Nextflow workflow that can be run as it’s own individual pipeline or included as part of another. 

CNViz is obvesously still underdevelopement, it's current version can be found as part of [this pipeline](https://github.com/ryanlayerlab/layer_lab_chco).

<img src="https://github.com/MSBradshaw/CNViz/blob/main/CNViz.png?raw=true" width="600"/>
