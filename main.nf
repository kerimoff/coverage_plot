#!/usr/bin/env nextflow
/*
========================================================================================
                          eQTL-Catalogue/qtlmap
========================================================================================
 eQTL-Catalogue/qtlmap Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/eQTL-Catalogue/qtlmap
----------------------------------------------------------------------------------------
*/
nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'

     eQTL-Catalogue/qtlmap v${workflow.manifest.version}
    =======================================================

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf\
     -profile tartu_hpc\
     --studyFile testdata/multi_test.tsv\
     --vcf_has_R2_field FALSE\
     --varid_rsid_map_file testdata/varid_rsid_map.tsv.gz\
     --n_batches 25

    Mandatory arguments:
      --studyFile                   Path to the TSV file containing pipeline inputs (VCF, expression matrix, metadata)

    Executions:
      --n_batches                   Number of parallel batches to run QTL Mapping per sample, must exceed the number of chromosomes (default: 400)
      --vcf_has_R2_field            Does the genotype VCF file contain R2 value in the INFO field? (default: true)
      --run_permutation             Calculate permuation p-values for each phenotype group (group_id in the phenotype metadata file) (default: false)
      --run_nominal                 Calculate nominal p-values for each phenotype group (group_id in the phenotype metadata file) (default: true)
      --n_permutations              Number of permutations to be performed per gene when run_permutation = true (default: 1000)

    QTL mapping:
      --cis_window                  The window where to search for associated variants around the phenotype (default: 1000000)
      --n_geno_pcs                  Number of genotype matrix principal components included as covariates in QTL analysis (default: 6).
      --n_pheno_pcs                 Number of phenotype matrix principal components included as covariates in QTL analysis (default: 6).
      --mincisvariant               Minimal numner of variants needed to be found in cis_window of each phenotype (default: 5)
      --covariates                  Comma-separated list of additional covariates included in the analysis (e.g. sex, age, batch). Columns with the exact same names should exist in the sample metadata file. 

    Fine mapping (SuSiE)
      --run_susie                   Perform eQTL fine mapping with SuSiE
      --vcf_genotype_field          Field in the VCF file that is used to construct the dosage matrix. Valid options are GT and DS (default: GT). 
      --susie_skip_full             Avoid writing full SuSiE output to disk. (default: false).

    Format results:
      --reformat_sumstats          Add rsid and median TPM columns to the nominal summary statistics files and perform additional formatting to make the files compatible with the eQTL Catalogue (default: true)
      --varid_rsid_map_file         TSV file mapping variant ids in CHR_POS_REF_ALT format to rsids from dbSNP.

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

/*
 * Create a channel for input files
 */ 

// quant_method	name_of_study	qtl_group	susie_purity_filtered	sample_meta	vcf_file	bigwig_path	usage_matrix_norm	exon_summ_stats	all_summ_stats
Channel.fromPath(params.studyFile)
    .ifEmpty { error "Cannot find studyFile file in: ${params.studyFile}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.name_of_study, row.quant_method, row.qtl_group, file(row.susie_purity_filtered), file(row.sample_meta), file(row.bigwig_path), file(row.usage_matrix_norm), file(row.exon_summ_stats), file("${row.exon_summ_stats}.tbi"), file(row.all_summ_stats), file("${row.all_summ_stats}.tbi"), file(row.pheno_meta), file(row.scaling_factors) ]}
    .branch {
        ge: it[1] == "ge"
        exon: it[1] == "exon"
        tx: it[1] == "tx"
        txrev: it[1] == "txrev"
        leafcutter: it[1] == "leafcutter"
    }
    .set { study_file_ch }

Channel.fromPath(params.studyFile)
  .ifEmpty { error "Cannot find studyFile file in: ${params.studyFile}" }
  .splitCsv(header: true, sep: '\t', strip: true)
  .map{row -> [ row.name_of_study, file(row.vcf_file) ]}
  .distinct()
  .set { vcf_file_ch }

// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

eQTL-Catalogue/qtlmap v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']        = 'eQTL-Catalogue/recap_plot'
summary['Pipeline Version']     = workflow.manifest.version
summary['Run Name']             = custom_runName ?: workflow.runName
summary['Study file']           = params.studyFile

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "========================================="

include { recap_plot_ge } from './workflows/recap_plot_ge_wf'
include { recap_plot_tx } from './workflows/recap_plot_tx_wf'
include { recap_plot_txrev } from './workflows/recap_plot_txrev_wf'
include { recap_plot_exon } from './workflows/recap_plot_exon_wf'
include { recap_plot_leafcutter } from './workflows/recap_plot_leafcutter_wf'
include { tabix_index } from './modules/utils'

workflow {
    tabix_index(vcf_file_ch)
    recap_plot_ge(study_file_ch.ge.combine(tabix_index.out, by: 0).distinct().view())
    recap_plot_tx(study_file_ch.tx.combine(tabix_index.out, by: 0).distinct())
    recap_plot_txrev(study_file_ch.txrev.combine(tabix_index.out, by: 0).distinct())
    recap_plot_exon(study_file_ch.exon.combine(tabix_index.out, by: 0).distinct())
    recap_plot_leafcutter(study_file_ch.leafcutter.combine(tabix_index.out, by: 0).distinct())
}

workflow.onComplete {
    log.info "[eQTL-Catalogue/recap_plot] Pipeline Complete"
}
