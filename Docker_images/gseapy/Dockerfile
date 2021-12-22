# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
# This is based on my scanpy:1.4.6 container - just install gseapy on top 
ARG BASE_CONTAINER=evafast1/scanpy:1.4.6

FROM $BASE_CONTAINER

LABEL maintainer="Eva Fast <evaisfast@gmail.com>"

USER root

# install with pip
RUN pip install --no-cache-dir\
  gseapy 
  
USER root
