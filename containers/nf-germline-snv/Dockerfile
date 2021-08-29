# continuumio/miniconda3:4.10.3 linux/amd64
FROM continuumio/miniconda3@sha256:8ba38cba2d88e3e094e3e3c32f1daaacb57a331aaf17a5ee0194c30d774903e0

LABEL description="Dockerfile containing all the requirements for the nf-germline-snv pipeline" \
      author="dembra96@gmail.com"

ARG ENV_NAME="base"

USER root

COPY environment.yml /
RUN conda env update -n ${ENV_NAME} -f environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/${ENV_NAME}/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name ${ENV_NAME} > ${ENV_NAME}_exported.yml

# Initialise bash for conda
RUN conda init bash

WORKDIR /data/

