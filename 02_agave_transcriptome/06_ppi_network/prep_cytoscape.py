# Protein-protein interaction network analysis (Methods 5.2.6): prep cytoscape.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import os, pandas as pd
D = "[diretório_transcriptoma]"
OUT = f"{D}/cytoscape_pc"; os.makedirs(OUT, exist_ok=True)
CAT = {"c_estresse": "Resposta a estresse", "c_proteinas": "Controle de qualidade", "c_traducao": "Tradução",
       "c_carboidratos": "Carboidratos", "c_dna": "Replicação e reparo"}
def nos(df):
    df = df.copy(); df.index.name = "id"
    c = df["cats"].fillna("")
    for k, v in CAT.items(): df[k] = c.str.contains(v, regex=False).astype(int)
    df["c_outras"] = (df[list(CAT)].sum(axis=1) == 0).astype(int)
    keep = ["rotulo", "direction", "log2FoldChange", "module", "x", "y"] + [x for x in ("organela",) if x in df] + list(CAT) + ["c_outras"]
    return df[keep]
def arestas(e, a, b, w, n):
    e = e.rename(columns={a: "source", b: "target", w: "peso"})[["source", "target", "peso"]]
    d = n["direction"]; e["tipo"] = ["-".join(sorted([d[s], d[t]])) for s, t in zip(e.source, e.target)]
    return e
nA = nos(pd.read_csv(f"{D}/rede_visual_nos.tsv", sep="\t", index_col=0))
eA = arestas(pd.read_csv(f"{D}/rede_visual_arestas.tsv", sep="\t"), "ref_A", "ref_B", "score_permitido", nA)
nB = nos(pd.read_csv(f"{D}/rede_consenso_K3_nos.tsv", sep="\t", index_col=0))
eB = arestas(pd.read_csv(f"{D}/rede_consenso_K3_arestas.tsv", sep="\t"), "A", "B", "K", nB)
for nome, n, e in (("A_abacaxi_dendrobium", nA, eA), ("B_consenso_K3", nB, eB)):
    n.to_csv(f"{OUT}/{nome}_nos.tsv", sep="\t"); e.to_csv(f"{OUT}/{nome}_arestas.tsv", sep="\t", index=False)
    print(nome, len(n), "nós,", len(e), "arestas |", e.tipo.value_counts().to_dict())
