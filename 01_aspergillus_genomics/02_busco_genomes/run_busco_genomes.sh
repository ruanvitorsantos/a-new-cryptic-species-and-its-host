#!/usr/bin/env bash
# BUSCO on genome assemblies, lineage aspergillus_odb12, and extraction of
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

# single-copy BUSCO proteins (Section 4.2.2).
# Usage: run from a directory containing the genomes (*.fna).

# (once) download the lineage dataset
# busco --download aspergillus_odb12

for i in *.fna
do
    SAMPLE=$(echo ${i} | sed "s/\.fna//")
    echo ${SAMPLE}

    busco -i ${SAMPLE}.fna \
          -l [diretório_busco_aspergillus_odb12]/ \
          -o ${SAMPLE} -m genome -c 35 -f

    # pool the single-copy BUSCO proteins and prefix each header with the strain identifier
    cat ${SAMPLE}/run_aspergillus_odb12/busco_sequences/single_copy_busco_sequences/*faa > ${SAMPLE}.faa
    sed -i 's,>,>'"${SAMPLE}"'_,' "${SAMPLE}.faa"
done
