#!/usr/bin/env nextflow
nextflow.enable.dsl=2

Channel
    .fromPath(params.mane_transcript_gene_map, checkIfExists: true)
    .set { mane_transcript_gene_map_ch }

Channel
    .fromPath(params.mane_gtf_file, checkIfExists: true)
    .set { mane_gtf_file_ch }

Channel
    .fromPath(params.tx_gtf_file, checkIfExists: true)
    .set { tx_gtf_file_ch }

include { generate_recap_plot_tx } from  '../modules/generate_plots'
include { split_into_batches } from  '../modules/utils'

workflow recap_plot_tx {
    take: 
    study_tsv_inputs_ch
    
    main:
    split_into_batches(study_tsv_inputs_ch)

    generate_recap_plot_tx(
        split_into_batches.out.study_tsv_inputs_ch.combine(split_into_batches.out.susie_batches.transpose(), by:[0,1,2]),
        mane_transcript_gene_map_ch.collect(),
        mane_gtf_file_ch.collect(),
        tx_gtf_file_ch.collect()
    )
}
