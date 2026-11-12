# Protein-protein interaction network analysis (Methods 5.2.6): conexoes.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import os, numpy as np, pandas as pd
D = "[diretório_transcriptoma]"
O = f"{D}/rede_abacaxi"
TUDO = os.environ.get("TUDO", "0") == "1"
OK = ["neighborhood", "fusion", "cooccurence", "coexpression", "experimental"]
NOME = {"neighborhood": "vizinhança gênica", "fusion": "fusão gênica", "cooccurence": "co-ocorrência", "coexpression": "coexpressão", "experimental": "experimentos"}
n = pd.read_csv(f"{O}/nos_completa_simbolos.tsv" if os.path.exists(f"{O}/nos_completa_simbolos.tsv") else f"{O}/nos_completa.tsv", sep="\t", index_col=0)
if "simbolo" not in n: n["simbolo"] = n.rotulo
e = pd.read_csv(f"{O}/{'arestas_completa' if TUDO else 'arestas'}.tsv", sep="\t")
if not TUDO: n = n[n.na_figura] if "na_figura" in n else n
need = set(e.source) | set(e.target)
inf = pd.read_csv(f"{D}/string_db/4615.protein.info.v12.0.txt.gz", sep="\t", header=0, usecols=[0, 3], names=["id", "ann"]).set_index("id")["ann"]
ch = {}
for c in pd.read_csv(f"{D}/string_db/4615.protein.links.detailed.v12.0.txt.gz", sep=" ", chunksize=2_000_000):
    c = c[c.protein1.isin(need) & c.protein2.isin(need) & (c.protein1 < c.protein2)]
    for r in c.to_dict("records"):
        ch[(r["protein1"], r["protein2"])] = {k: r[k] / 1000 for k in OK}
srt = lambda a, b: (a, b) if a < b else (b, a)
rows = []
for r in e.to_dict("records"):
    a, b = r["source"], r["target"]
    if a not in n.index or b not in n.index: continue
    d = ch.get(srt(a, b), {k: 0 for k in OK})
    dom = max(d, key=d.get) if max(d.values()) > 0 else "?"
    na, nb = n.loc[a], n.loc[b]
    rows.append(dict(gene_A=na.simbolo, anotacao_A=str(inf.get(a, ""))[:90], dir_A=na.direction, log2FC_A=round(na.log2FoldChange, 2), modulo_A=na.module,
                     gene_B=nb.simbolo, anotacao_B=str(inf.get(b, ""))[:90], dir_B=nb.direction, log2FC_B=round(nb.log2FoldChange, 2), modulo_B=nb.module,
                     tipo=r["tipo"], score=round(r["peso"], 3), canal_dominante=NOME.get(dom, "?"),
                     **{NOME[k]: round(d[k], 2) for k in OK}, id_A=a, id_B=b, transcrito_A=na.t, transcrito_B=nb.t))
C = pd.DataFrame(rows).sort_values(["tipo", "score"], ascending=[False, False])
suf = "_rede_completa" if TUDO else ""
C.to_csv(f"{O}/conexoes{suf}.tsv", sep="\t", index=False)
print("Conexões:", len(C), "|", C.tipo.value_counts().to_dict())
print("Canal dominante:", C.canal_dominante.value_counts().to_dict())
g = []
for i, r in n.loc[n.index.isin(need)].iterrows():
    x = C[(C.id_A == i) | (C.id_B == i)]
    par = [(row.gene_B, row.dir_B) if row.id_A == i else (row.gene_A, row.dir_A) for row in x.itertuples()]
    g.append(dict(gene=r.simbolo, direcao=r.direction, log2FC=round(r.log2FoldChange, 2), modulo=r.module, organela=r.get("organela", ""), conexoes=len(par),
                  com_up=sum(d == "up" for _, d in par), com_down=sum(d == "down" for _, d in par),
                  parceiros="; ".join(f"{p}({d[0]})" for p, d in sorted(par, key=lambda z: z[0])), id=i, transcrito=r.t))
G = pd.DataFrame(g).sort_values("conexoes", ascending=False)
G.to_csv(f"{O}/genes_resumo{suf}.tsv", sep="\t", index=False)
print("Genes up com conexão up-up:", int(((G.direcao == "up") & (G.com_up > 0)).sum()), "de", int((G.direcao == "up").sum()))
print(G.head(12)[["gene", "direcao", "conexoes", "com_up", "com_down"]].to_string(index=False))
print("Gravados:", f"conexoes{suf}.tsv e genes_resumo{suf}.tsv")
