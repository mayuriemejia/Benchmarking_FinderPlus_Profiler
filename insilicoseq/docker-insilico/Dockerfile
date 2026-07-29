FROM python:3.10-slim

# 1. Instalamos procps para que Nextflow pueda ejecutar 'ps' y medir métricas
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*

# 2. Instalamos InSilicoSeq
RUN pip install --no-cache-dir insilicoseq==2.0.1

# Dejamos que Nextflow controle el punto de entrada (sin ENTRYPOINT estricto)
