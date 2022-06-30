FROM bioconductor/bioconductor_docker:RELEASE_3_13
LABEL authors="Kaur Alasoo" \
      description="Docker image containing all requirements for susieR fine mapping"

RUN R -e "BiocManager::install(c('dplyr', 'optparse', 'readr', 'GenomicRanges', 'seqminer', 'rlist', 'purrr', 'rtracklayer', 'assertthat', 'cowplot'))"
RUN R -e "devtools::install_github(repo = 'kauralasoo/wiggleplotr@v22.06.1')"



