# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
# Customized using Jupyter Notebook R Stack https://github.com/jupyter/docker-stacks/tree/master/r-notebook
# Instructions for how to install from this Github source https://github.com/Bioconductor-notebooks/Identification-of-Differentially-Expressed-Genes-for-Ectopic-Pregnancy/blob/master/Dockerfile

ARG BASE_CONTAINER=jupyter/r-notebook
FROM $BASE_CONTAINER

MAINTAINER Eva Fast <evaisfast@gmail.com>

USER root

# R pre-requisites
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    fonts-dejavu \
    gfortran \
    pkg-config \
    zlib1g \
    gcc && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER $NB_USER

# R packages

# install in conda (bioconda channel)
RUN conda config --add channels bioconda && \
    conda config --add channels r && \
    conda config --add channels conda-forge && \
    conda config --add channels anaconda && \
    conda install -y \
    r-seurat\
    r-biocmanager \
    r-hdf5r \
    bioconductor-genomeinfodb \
    bioconductor-biovizbase \
    bioconductor-qvalue \
    bioconductor-ggbio \
    bioconductor-jaspar2018 \
    bioconductor-motifmatchr \
    bioconductor-chromvar \
    bioconductor-ensdb.mmusculus.v79 \
    bioconductor-go.db \
    bioconductor-tfbstools \
    bioconductor-bsgenome.mmusculus.ucsc.mm10 \
    bioconductor-diffloop \
    bioconductor-singlecellexperiment \
    cairo \
    pkgconfig \
    htslib \
    && \
    conda clean -tipy


# install packages from cran

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "remove.packages('data.table')"
RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "install.packages('rlist')"
RUN Rscript -e "install.packages('svglite')"
RUN Rscript -e "install.packages('devtools')"
RUN R -e "devtools::install_github('Rdatatable/data.table')"
RUN R -e "devtools::install_github('timoast/signac@v0.2.5')"


WORKDIR /home/jovyan
ADD . /home/jovyan
