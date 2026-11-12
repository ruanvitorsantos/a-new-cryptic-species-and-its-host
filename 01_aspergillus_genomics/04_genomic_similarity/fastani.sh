#!/usr/bin/env bash
# Optional cross-check with FastANI (not reported in the manuscript).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


conda install -c bioconda fastani -y
ls *.fasta > lista_genomas.txt
fastANI --ql lista_genomas.txt --rl lista_genomas.txt -o resultados_fastANI.txt
