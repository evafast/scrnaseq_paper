# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
ARG BASE_CONTAINER=jupyter/scipy-notebook

FROM $BASE_CONTAINER

LABEL maintainer="Eva Fast <evaisfast@gmail.com>"

USER root

# install in conda (bioconda channel)
RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y \
    python-igraph \
    louvain \
    leidenalg \
    ipywidgets \
    mnnpy \
    && \
    conda clean -tipy

# install with pip
RUN pip install --no-cache-dir\
  scanpy \
  gprofiler-official \
  UpSetPlot \
  'cellbrowser==0.7.7' \
  openpyxl \
  harmonypy \
  scanorama \
  h5py

USER root
