#!/usr/bin/env nextflow
nextflow.enable.dsl=2

Channel
    .fromPath(params.mane_transcript_gene_map, checkIfExists: true)
    .set { mane_transcript_gene_map_ch }

Channel
    .fromPath(params.mane_gtf_file, checkIfExists: true)
    .set { mane_gtf_file_ch }

include { generate_recap_plot_ge } from  '../modules/generate_plots'

workflow recap_plot_ge {
    take: 
    study_tsv_inputs_ch
    
    main:
    generate_recap_plot_ge(
        study_tsv_inputs_ch,
        mane_transcript_gene_map_ch.collect(),
        mane_gtf_file_ch.collect()
    )
}
