#!/usr/bin/env bash
# Reference proteomes (STRING v12.0) and canonical proteins for the DIAMOND searches used in the network analysis.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

set -e
D=[diretório_transcriptoma]
cd $D; mkdir -p ref5; cd ref5
zcat ../string_db5/{3702,4577,39947,4558,214687}.protein.sequences.v12.0.fa.gz > proteomas5.fa
grep -c ">" proteomas5.fa
cut -f1 ../diamond_agave_canonical_degs_best_hits.tsv | sort -u > ids_canonicos.txt
python3 - << 'PY'
ids = set(open("ids_canonicos.txt").read().split()); w = False; n = 0
with open("canonicos.faa", "w") as o:
    for l in open("../protein/Agave_sisalana_transcriptoma_sem_rRNA.fasta.transdecoder.pep"):
        if l[0] == ">": w = l[1:].split()[0] in ids; n += w
        if w: o.write(l)
print("Sequências canônicas extraídas:", n, "de", len(ids))
PY
diamond makedb --in proteomas5.fa -d prot5 --threads 8
diamond blastp -d prot5 -q canonicos.faa -o hits5_all.tsv --evalue 1e-5 -k 100 --threads 8 --outfmt 6 qseqid sseqid evalue bitscore pident length
python3 - << 'PY'
import pandas as pd
h = pd.read_csv("hits5_all.tsv", sep="\t", header=None, names=["q", "s", "e", "b", "p", "l"])
h["tax"] = h.s.str.split(".").str[0]
nq = len(open("ids_canonicos.txt").read().split())
print(f"Canônicos com hit em qualquer espécie: {h.q.nunique()} de {nq}")
print("Com hit, por espécie:", h.groupby("tax").q.nunique().to_dict())
print("Melhor hit global, por espécie:", h.sort_values("b", ascending=False).drop_duplicates("q").tax.value_counts().to_dict())
PY
