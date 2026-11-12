# Protein-protein interaction network analysis (Methods 5.2.6): semorf todos.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import re, pandas as pd
D = "[diretório_transcriptoma]"
P = f"{D}/protein/Agave_sisalana_transcriptoma_sem_rRNA.fasta"
DEG = "[diretório_deseq2_results]/DEGs_todos_shrink.tsv"
final = set()
for l in open(f"{P}.transdecoder.pep"):
    if l[0] == ">":
        final.add(re.sub(r"\.p\d+$", "", l[1:].split()[0]))
dg = pd.read_csv(DEG, sep="\t").set_index("gene")
deg = dg[(dg.padj < 0.05) & (dg.log2FoldChange.abs() > 1)]
ids = set(deg.index[~deg.index.isin(final)])
w = False; n = 0
with open(f"{D}/ref7/semorf_todos.fa", "w") as o:
    for l in open(P):
        if l[0] == ">":
            w = l[1:].split()[0] in ids; n += int(w)
        if w:
            o.write(l)
up = int((deg.loc[list(ids), "log2FoldChange"] > 0).sum())
print("DEGs sem ORF:", len(ids), "| up:", up, "| down:", len(ids) - up, "| sequências escritas:", n)
