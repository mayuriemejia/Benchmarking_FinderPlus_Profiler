FROM python:3.11-slim

RUN pip install openpyxl pandas
RUN apt-get update -qq \
 && apt-get install -y -qq --no-install-recommends procps \
 && rm -rf /var/lib/apt/lists/*

