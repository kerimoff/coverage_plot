process generate_recap_plot_ge {
    tag "${name_of_study}_${qtl_group}_${quant_method}"
    label "process_medium"
    publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "${name_of_study}_${qtl_group}_${quant_method}/*"
    publishDir "${params.outdir}/logs/${name_of_study}_${qtl_group}_${quant_method}/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
    tuple val(name_of_study), val(quant_method), val(qtl_group), file(sample_meta), file(bigwig_path), file(usage_matrix_norm), file(exon_summ_stats), file(exon_summ_stats_index), file(all_summ_stats), file(all_summ_stats_index), file(phenotype_meta), file(scaling_factors), file(vcf_file), file(vcf_file_index), file(susie_purity_filtered)
    path mane_transcript_gene_map
    path mane_gtf_file

    output:
    path "${name_of_study}_${qtl_group}_${quant_method}/*"
    path "${name_of_study}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log"

    script:
    individual_boxplots_check = params.individual_boxplots ? "--individual_boxplots" : ""
    """
    Rscript $projectDir/bin/plot_cov_ge.R \
        --name_of_study $name_of_study \
        --qtl_group $qtl_group \
        --finemap_susie $susie_purity_filtered \
        --sample_meta $sample_meta \
        --phenotype_meta $phenotype_meta \
        --vcf_file $vcf_file \
        --bigwig_path $bigwig_path \
        --mane_transcript_gene_map $mane_transcript_gene_map \
        --gtf_file $mane_gtf_file \
        --div_scaling_factors $scaling_factors \
        --usage_matrix_norm $usage_matrix_norm \
        --exon_summ_stats $exon_summ_stats \
        --all_summ_stats $all_summ_stats \
        --outdir ${name_of_study}_${qtl_group}_${quant_method} \
        $individual_boxplots_check 

    cp .command.log ${name_of_study}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log
    """
}

process generate_recap_plot_tx {
    tag "${name_of_study}_${qtl_group}_${quant_method}"
    label "process_medium"
    publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "${name_of_study}_${qtl_group}_${quant_method}/*"
    publishDir "${params.outdir}/logs/${name_of_study}_${qtl_group}_${quant_method}/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
    tuple val(name_of_study), val(quant_method), val(qtl_group), file(sample_meta), file(bigwig_path), file(usage_matrix_norm), file(exon_summ_stats), file(exon_summ_stats_index), file(all_summ_stats), file(all_summ_stats_index), file(phenotype_meta), file(scaling_factors), file(vcf_file), file(vcf_file_index), file(susie_purity_filtered)
    path mane_transcript_gene_map
    path mane_gtf_file
    path tx_gtf_file

    output:
    path "${name_of_study}_${qtl_group}_${quant_method}/*"
    path "${name_of_study}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log"

    script:
    individual_boxplots_check = params.individual_boxplots ? "--individual_boxplots" : ""
    """
    Rscript $projectDir/bin/plot_cov_tx.R \
        --name_of_study $name_of_study \
        --qtl_group $qtl_group \
        --finemap_susie $susie_purity_filtered \
        --sample_meta $sample_meta \
        --phenotype_meta $phenotype_meta \
        --vcf_file $vcf_file \
        --bigwig_path $bigwig_path \
        --mane_transcript_gene_map $mane_transcript_gene_map \
        --gtf_file $mane_gtf_file \
        --div_scaling_factors $scaling_factors \
        --tx_gtf_file $tx_gtf_file \
        --usage_matrix_norm $usage_matrix_norm \
        --exon_summ_stats $exon_summ_stats \
        --all_summ_stats $all_summ_stats \
        --outdir ${name_of_study}_${qtl_group}_${quant_method} \
        $individual_boxplots_check 

    cp .command.log ${name_of_study}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log
    """
}

process generate_recap_plot_txrev {
    tag "${name_of_study}_${qtl_group}_${quant_method}"
    label "process_medium"
    publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "${name_of_study}_${qtl_group}_${quant_method}/*"
    publishDir "${params.outdir}/logs/${name_of_study}_${qtl_group}_${quant_method}/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
    tuple val(name_of_study), val(quant_method), val(qtl_group), file(sample_meta), file(bigwig_path), file(usage_matrix_norm), file(exon_summ_stats), file(exon_summ_stats_index), file(all_summ_stats), file(all_summ_stats_index), file(phenotype_meta), file(scaling_factors), file(vcf_file), file(vcf_file_index), file(susie_purity_filtered)
    path mane_transcript_gene_map
    path mane_gtf_file
    path txrev_gtf_file

    output:
    path "${name_of_study}_${qtl_group}_${quant_method}/*"
    path "${name_of_study}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log"

    script:
    individual_boxplots_check = params.individual_boxplots ? "--individual_boxplots" : ""
    """
    Rscript $projectDir/bin/plot_cov_txrev.R \
        --name_of_study $name_of_study \
        --qtl_group $qtl_group \
        --finemap_susie $susie_purity_filtered \
        --sample_meta $sample_meta \
        --phenotype_meta $phenotype_meta \
        --vcf_file $vcf_file \
        --bigwig_path $bigwig_path \
        --mane_transcript_gene_map $mane_transcript_gene_map \
        --gtf_file $mane_gtf_file \
        --div_scaling_factors $scaling_factors \
        --txrev_gtf_file $txrev_gtf_file \
        --usage_matrix_norm $usage_matrix_norm \
        --exon_summ_stats $exon_summ_stats \
        --all_summ_stats $all_summ_stats \
        --outdir ${name_of_study}_${qtl_group}_${quant_method} \
        $individual_boxplots_check 

        cp .command.log ${name_of_study}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log
    """
}

process generate_recap_plot_exon {
    tag "${name_of_study}_${qtl_group}_${quant_method}"
    label "process_medium"
    publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "${name_of_study}_${qtl_group}_${quant_method}/*"
    publishDir "${params.outdir}/logs/${name_of_study}_${qtl_group}_${quant_method}/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
    tuple val(name_of_study), val(quant_method), val(qtl_group), file(sample_meta), file(bigwig_path), file(usage_matrix_norm), file(exon_summ_stats), file(exon_summ_stats_index), file(all_summ_stats), file(all_summ_stats_index), file(phenotype_meta), file(scaling_factors), file(vcf_file), file(vcf_file_index), file(susie_purity_filtered)
    path mane_transcript_gene_map
    path mane_gtf_file

    output:
    path "${name_of_study}_${qtl_group}_${quant_method}/*"
    path "${name_of_study}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log"

    script:
    individual_boxplots_check = params.individual_boxplots ? "--individual_boxplots" : ""
    """
    Rscript $projectDir/bin/plot_cov_exon.R \
        --name_of_study $name_of_study \
        --qtl_group $qtl_group \
        --finemap_susie $susie_purity_filtered \
        --sample_meta $sample_meta \
        --phenotype_meta $phenotype_meta \
        --vcf_file $vcf_file \
        --bigwig_path $bigwig_path \
        --mane_transcript_gene_map $mane_transcript_gene_map \
        --gtf_file $mane_gtf_file \
        --div_scaling_factors $scaling_factors \
        --usage_matrix_norm $usage_matrix_norm \
        --exon_summ_stats $exon_summ_stats \
        --all_summ_stats $all_summ_stats \
        --outdir ${name_of_study}_${qtl_group}_${quant_method} \
        $individual_boxplots_check 

    cp .command.log ${name_of_study}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log
    """
}

process generate_recap_plot_leafcutter {
    tag "${name_of_study}_${qtl_group}_${quant_method}"
    label "process_medium"
    publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "${name_of_study}_${qtl_group}_${quant_method}/*"
    publishDir "${params.outdir}/logs/${name_of_study}_${qtl_group}_${quant_method}/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
    tuple val(name_of_study), val(quant_method), val(qtl_group), file(sample_meta), file(bigwig_path), file(usage_matrix_norm), file(exon_summ_stats), file(exon_summ_stats_index), file(all_summ_stats), file(all_summ_stats_index), file(phenotype_meta), file(scaling_factors), file(vcf_file), file(vcf_file_index), file(susie_purity_filtered)
    path mane_transcript_gene_map
    path mane_gtf_file

    output:
    path "${name_of_study}_${qtl_group}_${quant_method}/*"
    path "${name_of_study}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log"

    script:
    individual_boxplots_check = params.individual_boxplots ? "--individual_boxplots" : ""
    """
    Rscript $projectDir/bin/plot_cov_leafcutter.R \
        --name_of_study $name_of_study \
        --qtl_group $qtl_group \
        --finemap_susie $susie_purity_filtered \
        --sample_meta $sample_meta \
        --phenotype_meta $phenotype_meta \
        --vcf_file $vcf_file \
        --bigwig_path $bigwig_path \
        --mane_transcript_gene_map $mane_transcript_gene_map \
        --gtf_file $mane_gtf_file \
        --div_scaling_factors $scaling_factors \
        --usage_matrix_norm $usage_matrix_norm \
        --exon_summ_stats $exon_summ_stats \
        --all_summ_stats $all_summ_stats \
        --outdir ${name_of_study}_${qtl_group}_${quant_method} \
        $individual_boxplots_check 

    cp .command.log ${name_of_study}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log
    """
}
