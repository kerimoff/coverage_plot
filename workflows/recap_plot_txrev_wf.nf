#!/usr/bin/env nextflow
nextflow.enable.dsl=2

Channel
    .fromPath(params.txrev_pheno_meta_path, checkIfExists: true)
    .set { pheno_metadata_ch }

Channel
    .fromPath(params.mane_transcript_gene_map, checkIfExists: true)
    .set { mane_transcript_gene_map_ch }

Channel
    .fromPath(params.mane_gtf_file, checkIfExists: true)
    .set { mane_gtf_file_ch }

Channel
    .fromPath(params.txrev_gtf_file, checkIfExists: true)
    .set { txrev_gtf_file_ch }

include { generate_recap_plot_txrev } from  '../modules/generate_plots'

workflow recap_plot_txrev {
    take: 
    study_tsv_inputs_ch
    
    main:
    generate_recap_plot_txrev(
        study_tsv_inputs_ch,
        pheno_metadata_ch.collect(),
        mane_transcript_gene_map_ch.collect(),
        mane_gtf_file_ch.collect(),
        txrev_gtf_file_ch.collect()
    )
}
