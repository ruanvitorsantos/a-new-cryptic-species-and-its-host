#!/usr/bin/env bash
# Retrieval of reference marker sequences (e.g. calmodulin) and search against the genomes.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


# one GenBank accession per line
cat << 'EOF' > accessions.txt
EF661154
...
EOF

# download all markers in one FASTA
efetch -db nucleotide -id "$(paste -sd',' accessions.txt)" -format fasta > marcadores_nigri.fasta

# nucleotide-vs-nucleotide-translated search against a genome database
makeblastdb -in genomes.fna -dbtype nucl -out db_genomes -parse_seqids
tblastx -query marcadores_nigri.fasta -db db_genomes -out resultado_tblastx.tsv \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

# protein-level alternative with DIAMOND (database built beforehand from the marker proteins)
diamond blastp \
    -q <(cat arquivo.faa arquivo2.faa) \
    -d [diretório_marcadores_cam]/db.dmnd \
    -f 6 -p 4 --max-target-seqs 5 \
    -o resultados_cam.tsv
