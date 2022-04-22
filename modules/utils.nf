process tabix_index {
    tag "${name_of_study}"
    storeDir "${projectDir}/vcfs"

    input:
    tuple val(name_of_study), file(vcf_file)

    output:
    tuple val(name_of_study), file(vcf_file), file("${vcf_file}.tbi")

    script:
    """
    tabix $vcf_file
    """
}