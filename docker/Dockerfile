# Reproducible environment for ProxLabToolkit (R 4.4 / Bioconductor 3.20).
# Build from the repository root:  docker build -f docker/Dockerfile -t proxlabtoolkit:1.0.0 .
# Run:                             docker run --rm -v "$PWD":/work -w /work proxlabtoolkit:1.0.0 Rscript Examples/example_NegKTSwt_vs_Mock.R
FROM bioconductor/bioconductor_docker:RELEASE_3_20
RUN Rscript -e 'install.packages(c("dplyr","tidyr","purrr","tibble","data.table","readr","metap","imputeLCMD"), repos="https://cloud.r-project.org")' \
 && Rscript -e 'BiocManager::install(c("limma","impute","DEP","DEqMS","MSstats","proDA","ROTS","qvalue","MsCoreUtils"), ask=FALSE, update=FALSE)'
COPY . /opt/ProxLabToolkit
RUN R CMD INSTALL /opt/ProxLabToolkit && Rscript /opt/ProxLabToolkit/tests/test_core.R
WORKDIR /work
