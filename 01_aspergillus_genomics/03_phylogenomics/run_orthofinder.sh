#!/usr/bin/env bash
# OrthoFinder on the pooled single-copy BUSCO proteins (Section 4.2.3).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

# Input: one .faa per genome (headers prefixed with the strain ID) in [diretório_filogenia].

micromamba activate qc_genomas      # env with OrthoFinder, DIAMOND and MAFFT
cd [diretório_filogenia]

# -M msa: multiple-sequence-alignment gene trees / concatenated supermatrix (SpeciesTreeAlignment.fa)
# -A mafft: aligner;  -S diamond: all-vs-all search
orthofinder -f . -t 8 -a 4 -M msa -A mafft -S diamond
