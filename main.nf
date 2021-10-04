#!/usr/bin/env nextflow

nextflow.preview.dsl=2

params.name = 'CNViz'
params.tag = 'latest' // Default tag is latest, to be overwritten by --tag <version>

build_db = params.function == "build_db"? true : false

// include process and workflows from external files
include {extractSampleInfo; singleColumnFileToChannel} from './lib/cnviz_utilities'
include {wf_cnv_build_panel_db;wf_cnviz} from './lib/wf_cnviz'

// Check if genome exists in the config file
if (params.genomes && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the genomes.config file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// load the sample list
sample_list_path = params.cnviz_input
sample_list = file(sample_list_path)
ch_input_sample = extractSampleInfo(sample_list)

workflow{
    // create channel for the probes or other regions on interest
    ch_target_bed = Channel.value(file(params.target_bed))
    // load all cnv calls
    ch_all_cnv_calls = singleColumnFileToChannel(file(params.all_cnv_calls))
    if (build_db) {
        // create a reference panel
        wf_cnv_build_panel_db(ch_input_sample, ch_target_bed, ch_all_cnv_calls)
    }else{
        if(params.cnviz_ref_panel_db){
            // load all of the necessary files for CNViz that should be parameters
            ch_genes_file = Channel.value(file(params.genes_file))
            ch_cnviz_ref_panel_db = Channel.value(file(params.cnviz_ref_panel_db))
            gnomad_sv_file = Channel.value(file(params.gnomad_sv_file))
            gnomad_sv_file_tbi = Channel.value(file(params.gnomad_sv_file_tbi))
            gnomad_sv = Channel.value([gnomad_sv_file, gnomad_sv_file_tbi])
            // run CNViz for each proband
            wf_cnviz(ch_input_sample,ch_target_bed,ch_all_cnv_calls,ch_cnviz_ref_panel_db,ch_genes_file,gnomad_sv_file,gnomad_sv_file_tbi)
        }else{
            println('Missing reference panel: the option --cnviz_ref_panel_db is required if not using --function build_db')
        }
    }
}
