#!/usr/bin/env bash
# Alignment-free ANI approximation with Mash (Section 4.2.4; Supplementary Table 3).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

# k-mer size 21, sketch size 10,000; ANI ~ (1 - D) x 100.

micromamba activate mash_env        # micromamba create -n mash_env -c bioconda mash

# sketch of all genomes
mash sketch -k 21 -s 10000 -o sketches *.fasta

# all-vs-all distances
mash dist sketches.msh sketches.msh > mash_dist.tab

# distance -> approximate ANI (%)
awk 'BEGIN{OFS="\t"} {print $1, $2, (1-$3)*100}' mash_dist.tab > mash_ani.tab
