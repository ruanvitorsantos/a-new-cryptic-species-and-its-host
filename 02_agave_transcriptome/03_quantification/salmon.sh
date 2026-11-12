#!/usr/bin/env bash
# Transcript quantification with Salmon (Section 5.2.2).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


conda create -n salmon_1_10 -c conda-forge -c bioconda salmon=1.10.3 -y
conda activate salmon_1_10

# index of the Trinity assembly
salmon index -t trinity_plantae.fasta -i salmon_index_trinity

# one quantification per library (output: one folder with quant.sf per sample)
for R1 in [diretório_salmon]/sisal/*_matched_R1.fastq.gz; do
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
    BASENAME=$(basename "$R1" | sed 's/_R1\.fastq\.gz//')

    echo "Quantificando a amostra: $BASENAME"

    salmon quant -i salmon_index_trinity \
                 -l A \
                 -1 "$R1" \
                 -2 "$R2" \
                 --validateMappings \
                 -p 8 \
                 --writeMappings=sam \
                 -o "[diretório_salmon_quants]/${BASENAME}"
done
