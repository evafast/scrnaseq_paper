# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
ARG BASE_CONTAINER=jupyter/r-notebook
FROM $BASE_CONTAINER

LABEL maintainer="Eva Fast <evaisfast@gmail.com>"

USER $NB_UID

# install in conda (bioconda channel)
RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y \
    bioconductor-scran \
    bioconductor-mast \
    python-igraph \
    louvain \
    ipywidgets \
    'scipy==1.3.0' \
    'statsmodels==0.10.0*' \
    && \
    conda clean -tipy

# install with pip
RUN pip install --no-cache-dir\
  scanpy \
  leidenalg \
  gprofiler-official \
  anndata2ri \
  UpSetPlot \
  h5py

# in actual environment enable: "! jupyter nbextension enable --py widgetsnbextension"

USER $NB_UID
