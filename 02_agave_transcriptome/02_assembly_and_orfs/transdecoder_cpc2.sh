#!/usr/bin/env bash
# Coding-region prediction (TransDecoder) and coding-potential check (CPC2).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


conda install -c bioconda transdecoder -y

# 1) long ORFs;  2) ORF selection by log-likelihood score
TransDecoder.LongOrfs -t arquivo.fasta
TransDecoder.Predict  -t arquivo.fasta
# main outputs (written to the current directory):
#   arquivo.fasta.transdecoder.pep   (proteins)  <- input for InterProScan / BLAST / DIAMOND
#   arquivo.fasta.transdecoder.cds   (coding sequences)
#   arquivo.fasta.transdecoder.gff3  (coordinates)

# optional: BLASTp of the proteins against UniProt/Swiss-Prot
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# CPC2 (coding-potential assessment)
wget https://github.com/gao-lab/CPC2_standalone/archive/refs/heads/master.zip
unzip master.zip
cd CPC2_standalone-master
conda install biopython -y
python bin/CPC2.py -i seu_transcriptoma.fa -o resultados.txt
