workflow wf_cnv_build_proband_db{
    take: _bams_vcfs
    take: _probes
    take: _calls
    main: 
    _bams = _bams_vcfs.map{ids, vcf, tbi, bam, bai  -> [ids, ids, bam, bai]}   
    _vcfs = _bams_vcfs.map{ids,vcf, tbi, bam, bai  -> [ids, ids, vcf, tbi]}
    _id = _bams_vcfs.map{ids,vcf, tbi, bam, bai  -> ids}


    mosdepth(_bams)
    gzip_probes(_probes)
    count_reads(_bams)
    get_total_reads(count_reads.out.map{idP, idS, counts -> [counts]}.collect())
    get_probe_reads_per_million(mosdepth.out,get_total_reads.out,gzip_probes.out)

    _rpm = get_probe_reads_per_million.out.map{pid,sid,rpm,gz -> [sid,rpm,gz]}
    _bams_vcfs_rpms = _bams_vcfs.join(_rpm,remainder: true)

    _id = _bams_vcfs.map{ids,vcf, tbi, bam, bai  -> ids}
    
    compile_proband_db(_bams_vcfs_rpms,
                mosdepth.out,
                get_total_reads.out,
                _calls.collect())

    emit:
    mosdepth = mosdepth.out
    rpm = get_probe_reads_per_million.out
    read_counts = count_reads.out
    total_read_counts = get_total_reads.out
    db = compile_proband_db.out
}

process load_panel_db{
    label 'container_llab'
    label 'cpus_1'
    publishDir "${params.outdir}/CNV_Plotting/LoadTestingArea", mode: params.publish_dir_mode

    input:
        file(panel_dir)

    output:
        tuple file('ReferencePanel_RPM'),  file('ReferencePanel_Calls'),  file('ReferencePanel_BAM'),  file('ReferencePanel_VCF'), file('ReferencePanel_MosDepth')


    script:
    """
    cp -r ${panel_dir}/* ./
    """
}

workflow wf_cnv_build_panel_db{
    take: _bams_vcfs
    take: _probes
    take: _calls

    main:

    _bams = _bams_vcfs.map{ids, vcf, tbi, bam, bai  -> [ids, ids, bam, bai]}
    _vcfs = _bams_vcfs.map{ids,vcf, tbi, bam, bai  -> [ids, ids, vcf, tbi]}

    mosdepth(_bams)
    gzip_probes(_probes)
    count_reads(_bams)
    get_total_reads(count_reads.out.map{idP, idS, counts -> [counts]}.collect())
    get_probe_reads_per_million(mosdepth.out,get_total_reads.out,gzip_probes.out)
    compile_db("ReferencePanel",mosdepth.out.collect(),get_probe_reads_per_million.out.collect(),count_reads.out.collect(),get_total_reads.out.collect(),_vcfs.collect(),_calls.collect())

    emit:
    mosdepth = mosdepth.out
    rpm = get_probe_reads_per_million.out
    read_counts = count_reads.out
    total_read_counts = get_total_reads.out   
    db = compile_db.out
}

workflow wf_CNViz_compile{
    take: proband_db
    take: ref_panel_db
    take: genes_file

    main:
    // combine all of the RPM results
    combine_rpm_files(ref_panel_db.collect(),proband_db)
    get_regions_zscores(combine_rpm_files.out)

    get_ref_RPM_files(ref_panel_db.collect())
    get_proband_rpm_files(proband_db)

    ref_rpm = get_ref_RPM_files.out
    pro_rpm = get_proband_rpm_files.out

    println(ref_rpm.view{"ref_rpm: $it"})
    println(pro_rpm.view{"pro_rpm: $it"})

    joint = pro_rpm.join(get_regions_zscores.out)
    
    println(joint.view{"joint: $it"})

    // get_adj_zscore(combine_rpm_files.out.map{id,files -> [files]}, get_regions_zscores.out)
    get_adj_zscore(ref_rpm, joint)
    get_probe_cover_mean_std(combine_rpm_files.out)
    merge_adj_scores(get_adj_zscore.out.collect())

    
    // allele balance portion
    //collect_allele_counts(_bams_unfiltered,_probes,_ref,_ref_fai,_ref_dict)
    //agg_allele_counts(collect_allele_counts.out, _probes)
    //merge_all_allele_counts(agg_allele_counts.out.map{ idP, idS, counts -> [counts]}.collect(), _probes)
    label_exons(genes_file)

    emit:
    rpm = combine_rpm_files.out
    adj_probe_scores = merge_adj_scores.out
    labeled_exons = label_exons.out
    probe_cover_mean_std = get_probe_cover_mean_std.out

    
}

process combine_rpm_files{
    label 'container_llab'
    label 'cpus_1'
    publishDir "${params.outdir}/CNV_Plotting/All_RPM_${sid}", mode: params.publish_dir_mode    

    input:
        tuple file(rpm), file(calls), file(bam), file(vcf), file(mosdepth)
        tuple sid, file(pro_rpm), file(pro_calls), file(pro_bam), file(pro_vcf), file(pro_mosdepth)

    output:
        tuple sid, file("*rpm_rate.bed*")
        

    script:
    """
    cp ${rpm}/* ./
    cp ${pro_rpm}/* ./ 
    """
}

process get_proband_rpm_files{
    label 'container_llab'
    label 'cpus_1'
    publishDir "${params.outdir}/CNV_Plotting/All_RPM_${sid}", mode: params.publish_dir_mode    

    input:
        tuple sid, file(pro_rpm), file(pro_calls), file(pro_bam), file(pro_vcf), file(pro_mosdepth)

    output:
        tuple sid, file("*rpm_rate.bed*")
        

    script:
    """
    cp ${pro_rpm}/* ./ 
    """
}

process get_ref_RPM_files{
    label 'container_llab'
    label 'cpus_1'
    publishDir "${params.outdir}/CNV_Plotting/Ref_RPM", mode: params.publish_dir_mode    

    input:
        tuple file(rpm), file(calls), file(bam), file(vcf), file(mosdepth)

    output:
        file("*rpm_rate.bed*")
        

    script:
    """
    cp ${rpm}/* ./
    """
}

process compile_proband_db{
    //label 'container_llab'
    label 'cpus_1'

    // tag {idSample}
    publishDir "${params.outdir}/CNV_Plotting/DB_${sid}", mode: params.publish_dir_mode
    input:
        tuple sid, file(vcf), file(tbi), file(bam), file(bai), file(rpm), file(rpm_tbi)
        file(all_mosdepth) // mosdepth
        file(all_total_read_count)// total read count
        file(all_calls)

    output:
        tuple sid, file("${sid}_RPM"), file("${sid}_Calls"), file("${sid}_BAM"), file("${sid}_VCF"), file("${sid}_MosDepth")

    script:
    """
    mkdir ${sid}_Calls
    mv $all_calls ${sid}_Calls
    
    mkdir ${sid}_VCF
    mv $vcf ${sid}_VCF
    mv $tbi ${sid}_VCF
    
    mkdir ${sid}_BAM
    mv $bam ${sid}_BAM
    mv $bai ${sid}_BAM
    
    mkdir ${sid}_RPM
    mv $rpm ${sid}_RPM
    mv $rpm_tbi ${sid}_RPM
    
    mkdir ${sid}_MosDepth
    mv $all_mosdepth ${sid}_MosDepth
    
    mkdir ${sid}_TotalReads
    mv $all_total_read_count ${sid}_TotalReads
    """
}

process compile_db{
    //label 'container_llab'
    label 'cpus_1'

    // tag {idSample}
    publishDir "${params.outdir}/CNV_Plotting/DB_${save_name}", mode: params.publish_dir_mode
    input:
        val(save_name)
        file(all_mosdepth) // mosdepth
        file(all_rpm)// rpm
        file(all_read_count)// read count
        file(all_total_read_count)// total read count
        file(all_vcfs)
        file(all_calls)

    output:
        tuple file("${save_name}_RPM"), file("${save_name}_Calls"), file("${save_name}_BAM"), file("${save_name}_VCF"), file("${save_name}_MosDepth")

    script:
    """
    mkdir ${save_name}_Calls
    mv $all_calls ${save_name}_Calls
    
    mkdir ${save_name}_VCF
    mv $all_vcfs ${save_name}_VCF
    
    mkdir ${save_name}_BAM
    
    mkdir ${save_name}_RPM
    mv $all_rpm ${save_name}_RPM
    
    mkdir ${save_name}_MosDepth
    mv $all_mosdepth ${save_name}_MosDepth
    
    mkdir ${save_name}_TotalReads
    mv $all_total_read_count ${save_name}_TotalReads
    """
}

workflow wf_cnv_data_prepossessing{
    take: _bams_unfiltered
    take: _probes
    take: _ref
    take: _ref_fai
    take: _ref_dict
    take: _genes_file

    main:
    mosdepth(_bams_unfiltered)
    gzip_probes(_probes)
    count_reads(_bams_unfiltered)
    get_total_reads(count_reads.out.map{idP, idS, counts -> [counts]}.collect())
    get_probe_reads_per_million(mosdepth.out,get_total_reads.out,gzip_probes.out)
    exon_coverage_rates(mosdepth.out, gzip_probes.out)
    
    get_regions_zscores(get_probe_reads_per_million.out.map{idP, idS, bed_gz, bed_gz_tbi -> [bed_gz, bed_gz_tbi]}.collect())
    get_adj_zscore(get_probe_reads_per_million.out, get_regions_zscores.out)
    get_probe_cover_mean_std(get_probe_reads_per_million.out.map{idP, idS, bed_gz, bed_gz_tbi -> [bed_gz, bed_gz_tbi]}.collect())
    merge_adj_scores(get_adj_zscore.out.map{idP, idS, bed_gz, bed_gz_tbi -> [bed_gz, bed_gz_tbi]}.collect())


    // allele balance portion
    collect_allele_counts(_bams_unfiltered,_probes,_ref,_ref_fai,_ref_dict)
    agg_allele_counts(collect_allele_counts.out, _probes)
    merge_all_allele_counts(agg_allele_counts.out.map{ idP, idS, counts -> [counts]}.collect(), _probes)
    label_exons(_genes_file)

    emit:
    adj_probe_scores = merge_adj_scores.out
    //allele_balance = merge_all_allele_counts.out
    //labeled_exons = label_exons.out
    //probe_cover_mean_std = get_probe_cover_mean_std.out
    
} //  end of wf_cnv_coverage_depth

process mosdepth {
    // mosdepth needs to be added to the llab container
    //label 'container_llab'
    label 'cpus_1'

    // tag {idSample}
    publishDir "${params.outdir}/CNV_Plotting/${idSample}/Mosdepth", mode: params.publish_dir_mode
    input:
        tuple idPatient, idSample, file(bam), file(bai)

    output:
        tuple idPatient, idSample, file("${idSample}.per-base.bed.gz"), file("${idSample}.per-base.bed.gz.tbi")


    // when: 'haplotypecaller' in tools

    script:
    """
    mosdepth $idSample $bam
    tabix -p bed ${idSample}.per-base.bed.gz
    touch del.txt
    """
}

// Get number of reads per bam
process count_reads {
    label 'container_llab'
    label 'cpus_1'
    publishDir "${params.outdir}/CNV_Plotting/${idSample}/ReadCounts", mode: params.publish_dir_mode 
    
    input:
        tuple idPatient, idSample, file(bam), file(bai)
 
    output:
        tuple idPatient, idSample, file("${idSample}.num_reads.txt")
    
    script:
    """
    samtools view -c -F 260 $bam > ${idSample}.num_reads.txt
    """

}


process get_probe_reads_per_million{
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/${idSample}/RPM", mode: params.publish_dir_mode

    input:
    tuple idPatient, idSample, file(per_base_bed_gz), file(per_base_bed_gz_tbi)
    file(total_reads) // out put from get_total_reads
    tuple file(probes_gz), file(probes_gz_tbi)

    output:
    tuple idPatient, idSample, file("${idSample}.probe.rpm_rate.bed.gz"), file("${idSample}.probe.rpm_rate.bed.gz.tbi")

    script:
    """
    get_rpm_rates.py \
    -m $per_base_bed_gz \
    --regions_file $probes_gz \
    --num_reads $total_reads \
    | bgzip -c > ${idSample}.probe.rpm_rate.bed.gz
    tabix -p bed ${idSample}.probe.rpm_rate.bed.gz
    """
}

process get_total_reads{
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/TotalReads", mode: params.publish_dir_mode

    input:
    file(num_read_files) // all the out puts from get_probe_reads_per_million

    output:
        file("total_read.txt")

    script:
    """
    cat *.num_reads.txt > total_read.txt
    """
}

process gzip_probes{
    label 'container_llab'
    label 'cpus_1'

    input:
    file(probes)
    
    output:
    tuple file("${probes}.sorted.gz"), file("${probes}.sorted.gz.tbi")
    
    """
    bedtools sort -i $probes > ${probes}.sorted
    bgzip ${probes}.sorted
    tabix ${probes}.sorted.gz -p bed
    """
}

process exon_coverage_rates {
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/${idSample}/ExonCoverageRates", mode: params.publish_dir_mode

    input:
    tuple idPatient, idSample, file(per_base_bed_gz), file(per_base_bed_gz_tbi) // mosdepth output
    tuple file(probes_gz), file(probes_gz_tbi) // gzip_probes output 

    output:
    tuple idPatient, idSample, file("${idSample}.probe.coverage_rate.bed.gz"), file("${idSample}.probe.coverage_rate.bed.gz.tbi")
        
    script:
    """
    get_coverage_rates.py \
    -m $per_base_bed_gz \
    --exons_file $probes_gz \
    | bgzip -c > ${idSample}.probe.coverage_rate.bed.gz
    tabix -p bed ${idSample}.probe.coverage_rate.bed.gz
    """
}

process get_regions_zscores {

    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/RegionZscores/${sid}", mode: params.publish_dir_mode

    input:
    tuple sid, file(all_rpm_files)// all the outputs from exon_coverage_rates
    // tuple sid, file(proband_rpm_gz), file(proband_rpm_gz_tbi) // all the outputs from exon_coverage_rates

    output:
    // removing the sid for now, it will be needed though
    //tuple file("probe.cover.mean.stdev.bed.gz"), file("probe.cover.mean.stdev.bed.gz.tbi"), file("probe.cover.mean.stdev.bed")
    tuple sid, file("probe.cover.mean.stdev.bed.gz"), file("probe.cover.mean.stdev.bed.gz.tbi"), file("probe.cover.mean.stdev.bed")

    script:
    """
    get_regions_zscores.py -r ./ | bgzip -c > probe.cover.mean.stdev.bed.gz
    cp probe.cover.mean.stdev.bed.gz copy_probe.cover.mean.stdev.bed.gz
    gunzip copy_probe.cover.mean.stdev.bed.gz 
    cp copy_probe.cover.mean.stdev.bed probe.cover.mean.stdev.bed
    tabix -p bed probe.cover.mean.stdev.bed.gz
    """   

}

process get_probe_cover_mean_std {

    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/ProbeCoverMeanSTD", mode: params.publish_dir_mode

    input:
    tuple sid, file(all_rpm_files)// all the outputs from exon_coverage_rates
    //tuple sid, file(probe_coverage_rate_bed_gz), file(probe_coverage_rate_bed_gz_tbi)// all the outputs from exon_coverage_rates

    output:
    tuple file("${sid}_probe.cover.mean.stdev.bed.gz"), file("${sid}_probe.cover.mean.stdev.bed.gz.tbi"), file("${sid}_probe.cover.mean.stdev.bed")

    script:
    """
    get_regions_zscores.py -r ./ | bgzip -c > ${sid}_probe.cover.mean.stdev.bed.gz
    cp ${sid}_probe.cover.mean.stdev.bed.gz ${sid}_copy_probe.cover.mean.stdev.bed.gz
    gunzip ${sid}_copy_probe.cover.mean.stdev.bed.gz 
    cp ${sid}_copy_probe.cover.mean.stdev.bed ${sid}_probe.cover.mean.stdev.bed
    tabix -p bed ${sid}_probe.cover.mean.stdev.bed.gz
    """

}

process get_adj_zscore{

    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/${sid}/AdjZscore", mode: params.publish_dir_mode

    input:
    tuple file(ref_rpm_files) // combine_rpm_files output
    tuple sid, file(rpm), file(rpm_tbi), file(coverage_rate_bed_gz), file(coverage_rate_bed_gz_tbi), file(coverage_rate_bed) //out out from get_regions_zscores.out

    output:
    tuple sid, file("${sid}.adj_z.bed.gz"), file("${sid}.adj_z.bed.gz.tbi")

    script:
    """
    get_coverage_zscores.py \
        -r ${sid}.probe.rpm_rate.bed.gz \
        -s probe.cover.mean.stdev.bed.gz \
        | bgzip -c > ${sid}.adj_z.bed.gz
    tabix -p bed ${sid}.adj_z.bed.gz
    """
}

process get_adj_zscore_og{

    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/${sid}/AdjZscore", mode: params.publish_dir_mode

    input:
    tuple file(all_rpm_files) // combine_rpm_files output
    tuple sid, file(coverage_rate_bed_gz), file(coverage_rate_bed_gz_tbi), file(coverage_rate_bed) //out out from get_regions_zscores.out

    output:
    tuple sid, file("${sid}.adj_z.bed.gz"), file("${sid}.adj_z.bed.gz.tbi")

    script:
    """
    get_coverage_zscores.py \
        -r ${sid}.probe.rpm_rate.bed.gz \
        -s $coverage_rate_bed_gz \
        | bgzip -c > ${sid}.adj_z.bed.gz
    tabix -p bed ${sid}.adj_z.bed.gz
    """
}

process merge_adj_scores{
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/MergeAdjZscore", mode: params.publish_dir_mode

    input:
    file(all_files) // collected output from get_adj_zscore
    // tuple sid, file(adj_z_bed_gz), file(adj_z_bed_gz_tbi) // collected output from get_adj_zscore

    output:
    tuple file("adj_scores.bed.gz"), file("adj_scores.bed.gz.tbi")

    script:
    """
    merge_samples_adj_scores.py -r ./ | bgzip -c > adj_scores.bed.gz
    tabix -p bed adj_scores.bed.gz
    """
}


// Allele Balance Calculations

process collect_allele_counts {
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/${idSample}/collect_allele_counts", mode: params.publish_dir_mode

    input:
    tuple idPatient, idSample, file(bam), file(bai)
    file(probes)
    file(ref)    
    file(fai)
    file(dict)

    output:
    tuple idPatient, idSample, file("${idSample}.allele_count.tsv")

    script:
    """
    gatk CollectAllelicCounts \
     -I $bam \
     -R $ref \
     -L $probes \
     -O ${idSample}.allele_count.tsv
    """
}

process agg_allele_counts{
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/${idSample}/agg_allele_counts", mode: params.publish_dir_mode

    input:
    tuple idPatient, idSample, file(collected_allele_counts) // output from collect_allele_counts
    file(probes)

    output:
    tuple idPatient, idSample, file("${idSample}.agg.allele_count.bed")

    script:
    """
    cat $collected_allele_counts | awk '{print \$1"\t"\$2"\t"\$2+1"\t"\$3"\t"\$4"\t"\$5"\t"\$6}' > ${idSample}.allele_count.bed
    bgzip ${idSample}.allele_count.bed
    tabix -p bed ${idSample}.allele_count.bed.gz
    agg_single_allele_counts.py $probes ${idSample}.allele_count.bed.gz ${idSample} > ${idSample}.agg.allele_count.bed    
    """
}

process merge_all_allele_counts {
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/merge_all_allele_counts", mode: params.publish_dir_mode

    input:
    file(agged_allele_counts) // output from agg_allele_counts
    file(probes)

    output:
    tuple file("ab.sorted.tsv.gz"), file("ab.sorted.tsv.gz.tbi"), file("ab.header.tsv")
    
    script:
    """
    cat *.agg.allele_count.bed > all_samples.aggregate_probe_allele_counts.txt
    create_all_sample_allele_count_bed.py all_samples.aggregate_probe_allele_counts.txt $probes> ab.bed
    head -1 ab.bed > ab.header.tsv
    tail -n +2 ab.bed > ab.tsv
    bedtools sort -i ab.tsv > ab.sorted.tsv
    bgzip ab.sorted.tsv
    tabix -p bed ab.sorted.tsv.gz 
    """
}

process label_exons{
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/label_exons", mode: params.publish_dir_mode

    input:
    file(exons_gz)

    output:
    tuple file("labeled_exons.bed.gz"), file("labeled_exons.bed.gz.tbi")

    script:
    """
    zcat $exons_gz | label_exon_number.py > labeled_exons.bed
    bgzip labeled_exons.bed
    tabix -p bed labeled_exons.bed.gz
    """
}

process combine_savvy_calls{
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/Savvy_call_params", mode: params.publish_dir_mode

    input:
        file(beds)

    output:
        file("savvy_call_plotting_params.txt")

    script:
    """
    for i in *.bed; do
        cat \$i | awk -v VAR=\$i '{ gsub(".coverageBinner","",\$10); print \$1":"\$2"-"\$3"\t"\$4"\t"\$10"\t"VAR"\t"\$1"."\$2"-"\$3"\t"\$4}' >> tmp_savvy_call_plotting_params.txt
    done

    for i in {1..23}; do
        grep "^\$i" tmp_savvy_call_plotting_params.txt >> savvy_call_plotting_params.txt || echo "Not found"
    done

    touch thing.txt
    """
}

process testy2{
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/TESTAREA2", mode: params.publish_dir_mode

    input:
        tuple thing1, thing2, thing3, thing4

    output:
        file("testy.txt")

    """
    touch testy.txt
    echo $thing1 >> testy.txt
    echo "\n_\n"
    echo $thing2 >> testy.txt
    echo "\n_\n"
    echo $thing3 >> testy.txt
    echo "\n_\n"
    echo $thing4 >> testy.txt
    """

}

process cnv_plotter_og {
    label 'container_py3_pandas'
    label 'cpus_16'

    publishDir "${params.outdir}/CNV_Plotting/${idSample}/Plots", mode: params.publish_dir_mode

    input:
    tuple file(ab_bed_gz), file(ab_bed_gz_tbi), file(ab_header) // output from merge_all_allele_counts
    tuple file(adj_scores_bed_gz), file(adj_scores_bed_gz_tbi)
    // tuple all_caller_type,  all_vcf_idPatient, all_vcf_idSample, file(all_vcf), file(all_vcf_tbi) // should contain all vcfs in the named like ${idSample}.vcf
    file(all_vcf)
    tuple file(labeled_exons_bed_gz), file(labeled_exons_bed_gz_tbi) // output from label_exons
    tuple file(probe_cover_mean_stdev_bed_gz), file(probe_cover_mean_stdev_bed_gz_tbi), file(probe_cover_mean_stdev_bed)
    tuple region, sv_type, idSample, filename, region_without_colon
    file(beds)

    output:
    file("${idSample}.${region_without_colon}.${sv_type}.png")

    script:
    """
    for i in *.bed; do
        if [ "\$i" != "probe.cover.mean.stdev.bed" ]
        then 
            bedtools sort -i \${i} > \${i}.sorted 2>> error.txt
            bgzip \${i}.sorted 2>> error.txt
            tabix -p bed \${i}.sorted.gz 2>> error.txt
            echo "\${i}.sorted.gz" >> multiple_savvy_calls.txt 2>> error.txt
        fi
    done   


    cnv_plotter.py --sample $idSample \
    --vcf ${idSample}.vcf.gz \
    -o ${idSample}.${region_without_colon}.${sv_type}.png \
    --scores $adj_scores_bed_gz \
    --exons $labeled_exons_bed_gz \
    --window 100000 \
    --height 7 \
    --width 5 \
    --region $region \
    --title "${idSample} ${region} ${sv_type}" \
    --alt_allele_counts $ab_bed_gz \
    --alt_allele_headers $ab_header \
    --label_exons \
    --depth $probe_cover_mean_stdev_bed \
    --all_calls multiple_savvy_calls.txt 2>> error.txt

    """

}

process get_max_number_of_calls {
    label 'container_py3_pandas'
    cpus = 8
    memory = '64 GB' 

    publishDir "${params.outdir}/CNV_Plotting/${sid}/MaxCalls", mode: params.publish_dir_mode

    input:
    //tuple file(adj_scores_bed_gz), file(adj_scores_bed_gz_tbi)
    //file(all_vcf)
    tuple file(labeled_exons_bed_gz), file(labeled_exons_bed_gz_tbi) // output from label_exons
    file(all_probe_cover_files)
    tuple region, sv_type, sid, filename, region_without_colon
    file(beds)

    output:
    file("${sid}.${region_without_colon}.${sv_type}.max_num_calls.txt")

    script:
    """
    for i in *.bed; do
        if [[ "\$i" != *"probe.cover.mean.stdev.bed" ]]
        then
            bedtools sort -i \${i} > \${i}.sorted 2>> error.txt
            bgzip \${i}.sorted 2>> error.txt
            tabix -p bed \${i}.sorted.gz 2>> error.txt
            echo "\${i}.sorted.gz" >> multiple_savvy_calls.txt 2>> error.txt
        fi
    done

    get_max_number_of_calls.py --sample $sid \
    --exons $labeled_exons_bed_gz \
    --window 100000 \
    --region $region \
    --depth ${sid}_probe.cover.mean.stdev.bed \
    --all_calls multiple_savvy_calls.txt > ${sid}.${region_without_colon}.${sv_type}.max_num_calls.txt
    """
}

process find_max_of_maxes {
    label 'container_py3_pandas'
    label 'cpus_16'

    publishDir "${params.outdir}/CNV_Plotting/MaxOfMaxes", mode: params.publish_dir_mode

    input:
        val all_maxes

    output:
        file('the_max_num_calls.txt')

    script:
    """
    echo $all_maxes > the_max_num_calls.txt
    #echo "${all_maxes}[*]" | sort -nr | head -1 > the_max_num_calls.txt 
    """
}

process cnv_plotter {
    label 'container_py3_pandas'
    cpus = 8
    memory = '64 GB'

    publishDir "${params.outdir}/CNV_Plotting/${sid}/Plots", mode: params.publish_dir_mode

    input:
    tuple file(adj_scores_bed_gz), file(adj_scores_bed_gz_tbi)
    // tuple all_caller_type,  all_vcf_idPatient, all_vcf_idSample, file(all_vcf), file(all_vcf_tbi) // should contain all vcfs in the named like ${idSample}.vcf
    file(all_vcf)
    tuple file(labeled_exons_bed_gz), file(labeled_exons_bed_gz_tbi) // output from label_exons
    file(all_probe_cover_files)
    tuple region, sv_type, sid, filename, region_without_colon
    file(beds)
    file(max_num_calls)
    file(gnomad_sv)
    file(gnomad_sv_tbi)

    output:
    file("${sid}.${region_without_colon}.${sv_type}.png")

    script:
    """
    for i in *.bed; do
        if [[ "\$i" != *"probe.cover.mean.stdev.bed" ]]
        then
            bedtools sort -i \${i} > \${i}.sorted 2>> error.txt
            bgzip \${i}.sorted 2>> error.txt
            tabix -p bed \${i}.sorted.gz 2>> error.txt
            echo "\${i}.sorted.gz" >> multiple_savvy_calls.txt 2>> error.txt
        fi
    done

    cnv_plotter.py --sample $sid \
    --vcf ${sid}.vcf.gz \
    -o ${sid}.${region_without_colon}.${sv_type}.png \
    --scores $adj_scores_bed_gz \
    --exons $labeled_exons_bed_gz \
    --window 100000 \
    --height 7 \
    --width 5 \
    --region $region \
    --title "${sid} ${region} ${sv_type}" \
    --label_exons \
    --depth ${sid}_probe.cover.mean.stdev.bed \
    --max_num_calls $max_num_calls \
    --gnomad_sv $gnomad_sv \
    --all_calls multiple_savvy_calls.txt 2>> error.txt

    """

}

process extract_sample_specific_calls{
    publishDir "${params.outdir}/CNV_Plotting/FilteredCalls", mode: params.publish_dir_mode
    input:
        file(all_calls)
        val(ids)

    output:
        file("filtered_calls.bed")

    script:
    """
    cat *.bed > temp_filtered_calls.bed
    echo "$ids" > ids.txt
    extract_sample_specific_calls.py temp_filtered_calls.bed  ids.txt filtered_calls.bed
    """
}












