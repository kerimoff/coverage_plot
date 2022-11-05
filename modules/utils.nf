process tabix_index {
    tag "${name_of_study}"
    storeDir "${projectDir}/vcfs"
    container "quay.io/biocontainers/tabix:1.11--hdfd78af_0"

    input:
    tuple val(name_of_study), file(vcf_file)

    output:
    tuple val(name_of_study), file(vcf_file), file("${vcf_file}.tbi")

    script:
    """
    tabix $vcf_file
    """
}

process split_into_batches {
    tag "${name_of_study}_${qtl_group}_${quant_method}"
    label "process_low"

    input:
    tuple val(name_of_study), val(quant_method), val(qtl_group), file(susie_purity_filtered), file(sample_meta), file(bigwig_path), file(usage_matrix_norm), file(exon_summ_stats), file(exon_summ_stats_index), file(all_summ_stats), file(all_summ_stats_index), file(phenotype_meta), file(scaling_factors), file(vcf_file), file(vcf_file_index)

    output:
    tuple val(name_of_study), val(quant_method), val(qtl_group), file(sample_meta), file(bigwig_path), file(usage_matrix_norm), file(exon_summ_stats), file(exon_summ_stats_index), file(all_summ_stats), file(all_summ_stats_index), file(phenotype_meta), file(scaling_factors), file(vcf_file), file(vcf_file_index), emit: study_tsv_inputs_ch
    tuple val(name_of_study), val(quant_method), val(qtl_group), file("${name_of_study}_${qtl_group}_${quant_method}/*"), emit: susie_batches

    script:
    debug_mode_check = params.debug_mode ? "--debug_mode" : ""
    """
    Rscript $projectDir/bin/prepare_batches.R \
        --qtl_group $qtl_group \
        --quant_method $quant_method \
        --chunk_size ${params.chunk_size}\
        --finemap_susie $susie_purity_filtered \
        --phenotype_meta $phenotype_meta \
        --outdir ${name_of_study}_${qtl_group}_${quant_method} \
        $debug_mode_check 
    """
}