#!/usr/bin/env bash
# De novo transcriptome assembly with Trinity (Section 5.2.2).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


conda create -n trinity_env -y
conda activate trinity_env
conda config --env --add channels defaults
conda config --env --add channels bioconda
conda config --env --add channels conda-forge
conda config --env --set channel_priority strict
conda install -c bioconda trinity=2.15.1 -y

Trinity --seqType fq --max_memory 63G \
        --left  $(ls *fwd.fq.gz | paste -sd ',') \
        --right $(ls *rev.fq.gz | paste -sd ',') \
        --CPU 16 --normalize_reads --normalize_max_read_cov 50 \
        --output trinity_sisal

# assembly metrics (N50, etc.)
TrinityStats.pl trinity_sisal.Trinity.fasta
