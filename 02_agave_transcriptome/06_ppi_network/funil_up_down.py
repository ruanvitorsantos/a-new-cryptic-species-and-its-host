# Protein-protein interaction network analysis (Methods 5.2.6): funil up down.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import re, numpy as np, pandas as pd

D = "[diretório_transcriptoma]"
DEG = "[diretório_deseq2_results]/DEGs_todos_shrink.tsv"
Q = "[diretório_salmon_quants]"
PEP = f"{D}/protein/Agave_sisalana_transcriptoma_sem_rRNA.fasta.transdecoder.pep"
SD = ["ASD2_concatenado_S28_S17_S6", "ASD3_concatenado_S29_S18_S7", "ASD5_concatenado_S30_S19_S8",
      "ASD6_concatenado_S31_S20_S9", "ASD7_concatenado_S32_S21_S10", "ASD8_concatenado_S33_S22_S11"]
SS = ["ASS2_concatenado_S23_S12_S1", "ASS5_concatenado_S24_S13_S2", "ASS7_concatenado_S25_S14_S3",
      "ASS8_concatenado_S26_S15_S4", "ASS9_concatenado_S27_S16_S5"]
tid = lambda s: s.str.replace(r"\.p\d+$", "", regex=True)

dg = pd.read_csv(DEG, sep="\t").set_index("gene")
deg = dg[(dg.padj < 0.05) & (dg.log2FoldChange.abs() > 1)].copy()
deg["dir"] = np.where(deg.log2FoldChange > 0, "up", "down")

orf = set()
with open(PEP) as f:
    for l in f:
        if l[0] == ">": orf.add(re.sub(r"\.p\d+$", "", l[1:].split()[0]))
hits = set(tid(pd.read_csv(f"{D}/diamond_agave_best_hits.tsv", sep="\t", header=None, usecols=[0])[0]))
canon = set(tid(pd.read_csv(f"{D}/diamond_agave_canonical_degs_best_hits.tsv", sep="\t", header=None, usecols=[0])[0]))
net = pd.read_csv(f"{D}/rede_canais_permitidos_thr0.7.tsv", sep="\t")
rede = set(tid(pd.concat([net.agave_A, net.agave_B])))
col = pd.read_csv(f"{D}/rede_colapsada_nos.tsv", sep="\t", index_col=0)
colap = set(col.transcript_id)

print("Lendo TPM...")
need = set(deg.index); tp = {}
for s in SD + SS:
    q = pd.read_csv(f"{Q}/clean_{s}_matched/quant.sf", sep="\t", usecols=["Name", "TPM"])
    tp[s] = q[q.Name.isin(need)].set_index("Name")["TPM"]
tp = pd.DataFrame(tp)
tpm_ok = set(tp.index[(tp[SD].mean(axis=1) >= 1) | (tp[SS].mean(axis=1) >= 1)])

est = [("DEGs (padj<0,05; |log2FC|>1)", set(deg.index)), ("+ ORF predito", orf), ("+ hit no DIAMOND", hits),
       ("+ TPM >= 1 em D ou S", tpm_ok), ("+ conjunto canônico", canon), ("+ aresta em 0,7 (canais permitidos)", rede),
       ("+ rede colapsada (comp. >= 3 nós)", colap)]
atual = None; ant = None; linhas = []
for nome, s in est:
    atual = set(deg.index) if atual is None else atual & s
    d = deg.loc[list(atual), "dir"]
    up, dn = int((d == "up").sum()), int((d == "down").sum())
    r = ("-", "-") if ant is None else (f"{100 * up / ant[0]:.1f}", f"{100 * dn / ant[1]:.1f}")
    linhas.append((nome, up + dn, up, dn, f"{100 * up / (up + dn):.1f}", r[0], r[1])); ant = (up, dn)
print(pd.DataFrame(linhas, columns=["etapa", "n", "up", "down", "% up", "retenção up (%)", "retenção down (%)"]).to_string(index=False))

gene = lambda i: re.sub(r"_i\d+$", "", i)
for nome, ids in (("DEGs", set(deg.index)), ("conjunto canônico", canon & set(deg.index))):
    g = pd.DataFrame({"g": [gene(i) for i in ids], "d": deg.loc[list(ids), "dir"].values}).drop_duplicates()
    print(f"Genes únicos em {nome}: up = {(g.d == 'up').sum()} | down = {(g.d == 'down').sum()}")
