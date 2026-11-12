# Protein-protein interaction network analysis (Methods 5.2.6): blastx semorf.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import re, numpy as np, pandas as pd
D = "[diretório_transcriptoma]"
P = f"{D}/protein/Agave_sisalana_transcriptoma_sem_rRNA.fasta"
DEG = "[diretório_deseq2_results]/DEGs_todos_shrink.tsv"
final = set()
for l in open(f"{P}.transdecoder.pep"):
    if l[0] == ">":
        final.add(re.sub(r"\.p\d+$", "", l[1:].split()[0]))
dg = pd.read_csv(DEG, sep="\t").set_index("gene")
deg = dg[(dg.padj < 0.05) & (dg.log2FoldChange.abs() > 1)].copy()
deg["dir"] = np.where(deg.log2FoldChange > 0, "up", "down")
sem = deg[~deg.index.isin(final)]
rng = np.random.default_rng(42)
up = sem[sem.dir == "up"]
dn = sem[sem.dir == "down"]
amostra = list(up.sample(min(3000, len(up)), random_state=1).index) + list(dn.sample(min(3000, len(dn)), random_state=1).index)
pd.Series(amostra).to_csv("ref7/semorf_ids.txt", index=False, header=False)
print("Sem ORF: up =", len(up), "| down =", len(dn), "| amostra:", len(amostra))
