#!/usr/bin/env bash
# Post-assembly purity check: slice the assembly into 10-kb windows and BLAST them against NCBI nt.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


conda install -c bioconda seqkit -y

# 10-kb non-overlapping windows
seqkit sliding -W 10000 -s 10000 A._welwitschiae_CCMB_674.fasta > A._welwitschiae_CCMB_674_10kb.fasta

blastn -query A._welwitschiae_CCMB_674_10kb.fasta \
       -db [diretório_banco_nt]/nt \
       -outfmt 6 -out resultado_blastn.tsv -num_threads 16

# unique subject accessions of the hits (column 2); upload to NCBI Batch Entrez
# (https://www.ncbi.nlm.nih.gov/sites/batchentrez, Database = Nucleotide) to retrieve their taxonomy
cut -f2 resultado_blastn.tsv | sort -u > acessos_ncbi.txt
