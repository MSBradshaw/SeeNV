#!/usr/bin/env nextflow

nextflow.preview.dsl=2

params.name = 'CNViz'
params.tag = 'latest' // Default tag is latest, to be overwritten by --tag <version>

build_db = params.function == "build_db"? true : false

// include process and workflows from external files
include {extractSampleInfo; singleColumnFileToChannel} from './lib/utilities' 
include {wf_cnv_build_proband_db;wf_cnv_build_panel_db;load_panel_db;wf_CNViz_compile;combine_savvy_calls;get_max_number_of_calls;find_max_of_maxes;cnv_plotter;extract_sample_specific_calls} from './lib/wf_cnviz'

// Check if genome exists in the config file
if (params.genomes && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the genomes.config file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// load the sample list
sample_list_path = params.input
sample_list = file(sample_list_path)
ch_input_sample = extractSampleInfo(sample_list)

workflow{
    // create channel for the probes or other regions on interest
    ch_target_bed = Channel.value(file(params.target_bed))
    // load all cnv calls
    ch_all_cnv_calls = singleColumnFileToChannel(file(params.all_cnv_calls))
    if (build_db) {
        println("make a bd")
        // create the proband_db
        wf_cnv_build_panel_db(ch_input_sample, ch_target_bed, ch_all_cnv_calls)
    }else{
        if(params.cnviz_ref_panel_db){
            // load all of the necessary files for CNViz that should be parameters
            ch_genes_file = Channel.value(file(params.genes_file))
            ch_cnviz_ref_panel_db = Channel.value(file(params.cnviz_ref_panel_db))
            gnomad_sv_file = Channel.value(file(params.gnomad_sv_file))
            gnomad_sv_file_tbi = Channel.value(file(params.gnomad_sv_file_tbi))
            gnomad_sv = Channel.value([gnomad_sv_file, gnomad_sv_file_tbi])


            println("Nice benny")
            wf_cnv_build_proband_db(ch_input_sample, ch_target_bed, ch_all_cnv_calls)
            load_panel_db(ch_cnviz_ref_panel_db)
            wf_CNViz_compile(wf_cnv_build_proband_db.out.db, load_panel_db.out, ch_genes_file)

            ids = ch_input_sample.map{id,vcf, tbi, bam, bai  -> id}            
            extract_sample_specific_calls(ch_all_cnv_calls.collect(), ids.collect())

            combine_savvy_calls(extract_sample_specific_calls.out)

            ch_savvy_calls = combine_savvy_calls.out.splitText(){it.split("\t")}.map{ x -> [x[0],x[1],x[2],x[3],x[4]] }


            get_max_number_of_calls(
                wf_CNViz_compile.out.labeled_exons.collect(),
                wf_CNViz_compile.out.probe_cover_mean_std.collect(),
                ch_savvy_calls,
                ch_all_cnv_calls.collect())
            find_max_of_maxes(get_max_number_of_calls.out.splitText().map{it -> it.trim()}.collect())
            cnv_plotter( wf_CNViz_compile.out.adj_probe_scores,
                ch_input_sample.map{ids, vcf, tbi, bam, bai  -> [vcf,tbi]}.collect(),
                wf_CNViz_compile.out.labeled_exons,
                wf_CNViz_compile.out.probe_cover_mean_std.collect(),
                ch_savvy_calls,
                ch_all_cnv_calls.collect(),
                find_max_of_maxes.out,
                gnomad_sv_file,
                gnomad_sv_file_tbi)
        }else{
            println('No reference DB listed')
        }
    }
}
