process plot_from_data {
    tag "${rds_or_tar_file.simpleName}"
    label "process_low"
    publishDir "${params.outdir}/", mode: 'copy'

    input:
    path rds_or_tar_file

    output:
    path "plots_from_data/*"

    script:
    individual_boxplots_check = params.individual_boxplots ? "--individual_boxplots" : ""
    if (rds_or_tar_file.getName().contains(".tar.gz")) {
        """
        Rscript $projectDir/bin/plot_from_data.R \
            --tar_file $rds_or_tar_file \
            --outdir plots_from_data \
            $individual_boxplots_check 
        """
    } else {

        """
        Rscript $projectDir/bin/plot_from_data.R \
            --rds_file $rds_or_tar_file \
            --outdir plots_from_data \
            $individual_boxplots_check 
        """
    }
}
