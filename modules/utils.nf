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