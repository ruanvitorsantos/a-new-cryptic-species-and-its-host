#!/usr/bin/env bash
# Maximum-likelihood species tree from the OrthoFinder supermatrix (Section 4.2.3).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

# ModelFinder Plus (-m MFP), 1,000 ultrafast bootstrap replicates (-bb 1000).

iqtree -s SpeciesTreeAlignment.fa -m MFP -bb 1000 -nt 12

# Same run, low priority on a shared server:
# nice -n 19 ionice -c3 taskset -c 0-11 iqtree -s SpeciesTreeAlignment.fa -m MFP -bb 1000 -nt 16 -mem 50G
