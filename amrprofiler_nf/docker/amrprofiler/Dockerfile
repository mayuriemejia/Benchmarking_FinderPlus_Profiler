# (•˕ •マ.ᐟ
# Imagen Docker para AMRProfiler v1.0.0
# Mayurie Mejía | Trainee/Rookie | Sequentia Biotech

FROM continuumio/miniconda3:23.10.0-1

LABEL maintainer="Mayurie Mejía"
LABEL version="1.0.0"

WORKDIR /amrprofiler

RUN apt-get update && \
    apt-get install -y --no-install-recommends unzip && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

COPY amrprofiler_db/amrprofiler-main/ ./

RUN conda env create -f /amrprofiler/environment.yml && \
    conda clean -afy

ENV PATH="/opt/conda/envs/amrprofiler/bin:$PATH"
ENV CONDA_DEFAULT_ENV=amrprofiler