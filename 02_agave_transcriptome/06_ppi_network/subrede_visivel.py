# Protein-protein interaction network analysis (Methods 5.2.6): subrede visivel.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import re, numpy as np, pandas as pd, networkx as nx
from functools import lru_cache
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Patch
from matplotlib.collections import LineCollection
from matplotlib.colors import TwoSlopeNorm
from matplotlib.cm import ScalarMappable

D = "[diretório_transcriptoma]"
IPS = "[diretório_interproscan]/resultados/ASS_novo_anotacao.tsv"
N = 150; SEED = 42

nodes = pd.read_csv(f"{D}/rede_colapsada_nos.tsv", sep="\t", index_col=0)
edges = pd.read_csv(f"{D}/rede_colapsada_arestas.tsv", sep="\t")
G = nx.from_pandas_edgelist(edges, "ref_A", "ref_B", edge_attr="score_permitido")
m = nodes.loc[list(G.nodes())].copy()

# seleção
m["escore"] = m[["degree", "betweenness"]].rank(pct=True).mean(axis=1)
forced = set(m.index[m.direction == "up"])
tam = m.module.value_counts(); grandes = tam[tam >= 10]
quota = (max(N - len(forced), 0) * grandes / grandes.sum()).round().astype(int).clip(lower=1)
sel = set(forced)
for mod, q in quota.items():
    sel |= set(m[(m.module == mod) & ~m.index.isin(forced)].sort_values("escore", ascending=False).head(q).index)
S = G.subgraph(sel).copy()
for c in [c for c in nx.connected_components(S) if len(c) < 3]: S.remove_nodes_from(c)
m = m.loc[list(S.nodes())]
print(f"Sub-rede: {S.number_of_nodes()} nós, {S.number_of_edges()} arestas | up que permaneceram: {(m.direction == 'up').sum()} de {len(forced)} | down: {(m.direction == 'down').sum()}")

# rótulos (anotação do STRING)
info = pd.concat([pd.read_csv(f"{D}/string_db/{t}.protein.info.v12.0.txt.gz", sep="\t", usecols=[0, 3], names=["id", "ann"], header=0) for t in (4615, 906689)]).drop_duplicates("id").set_index("id")["ann"]
def rot(r):
    a = info.get(r, "")
    a = re.split(r"[,;(]", a)[0].strip() if isinstance(a, str) else ""
    return a[:16] if a else "UNK"
m["rotulo"] = [rot(r) for r in m.index]

# GO -> categorias (raízes e descendentes)
par = {}; cur = None
for l in open(f"{D}/go-basic.obo"):
    l = l.strip()
    if l.startswith("["): cur = None
    if l == "[Term]": cur = "?"
    elif l.startswith("id: GO:") and cur == "?": cur = l[4:]; par[cur] = set()
    elif cur and cur != "?" and l.startswith("is_a: GO:"): par[cur].add(l[6:16])
    elif cur and cur != "?" and l.startswith("relationship: part_of GO:"): par[cur].add(l[22:32])
@lru_cache(None)
def anc(g): return frozenset({g}).union(*[anc(p) for p in par.get(g, ())]) if g in par else frozenset({g})
CAT = {"Resposta a estímulo e estresse": ["GO:0050896"],
       "Controle de qualidade, dobramento e degradação de proteínas": ["GO:0006457", "GO:0006508", "GO:0043161", "GO:0000502"],
       "Tradução e ribossomo": ["GO:0006412", "GO:0042254", "GO:0005840"],
       "Metabolismo de carboidratos e crescimento": ["GO:0005975", "GO:0040007"],
       "Replicação e reparo de DNA": ["GO:0006260", "GO:0006281"]}
COR = dict(zip(CAT, ["#D81B60", "#F2C230", "#8E5AB5", "#4CAF50", "#2B7BC9"]))
ips = pd.read_csv(IPS, sep="\t", header=None, usecols=[0, 13], names=["p", "go"], dtype=str)
ips = ips[ips.p.isin(set(m.node))]
gos = ips.groupby("p")["go"].apply(lambda s: set(re.findall(r"GO:\d{7}", " ".join(s.dropna()))))
def cats(p):
    A = set()
    for g in gos.get(p, set()): A |= anc(g)
    return [c for c, rs in CAT.items() if A & set(rs)]
m["categorias"] = [";".join(cats(p)) for p in m.node]
print("Nós por categoria:"); [print(f"  {c}: {int(m.categorias.str.contains(re.escape(c)).sum())}") for c in CAT]
print(f"Nós sem categoria: {int((m.categorias == '').sum())} | rótulo UNK: {int((m.rotulo == 'UNK').sum())}")

# tabelas para o Cytoscape
m.drop(columns=["abs"], errors="ignore").to_csv(f"{D}/subrede_{N}_nos.tsv", sep="\t")
nx.to_pandas_edgelist(S).rename(columns={"source": "ref_A", "target": "ref_B"}).to_csv(f"{D}/subrede_{N}_arestas.tsv", sep="\t", index=False)

# figura
pos = nx.spring_layout(S, seed=SEED, k=1.6 / np.sqrt(len(S)), iterations=300, weight="score_permitido")
xy = np.array(list(pos.values())); span = max(np.ptp(xy[:, 0]), np.ptp(xy[:, 1])); r = 0.022 * span
fig, ax = plt.subplots(figsize=(15, 15))
sg = [(pos[a], pos[b]) for a, b in S.edges()]
lw = [0.4 + 2.2 * (S[a][b]["score_permitido"] - 0.7) / 0.3 for a, b in S.edges()]
ax.add_collection(LineCollection(sg, colors="#8A8A8A", linewidths=lw, alpha=0.5, zorder=1))
norm = TwoSlopeNorm(vmin=-6, vcenter=0, vmax=6); cmap = plt.get_cmap("RdBu_r")
for n in S.nodes():
    x, y = pos[n]; row = m.loc[n]
    ax.add_patch(Circle((x, y), r * 0.8, fc=cmap(norm(np.clip(row.log2FoldChange, -6, 6))), ec="none", zorder=2))
    cs = [c for c in row.categorias.split(";") if c]
    if cs:
        for i, c in enumerate(cs):
            ax.add_patch(Wedge((x, y), r, 360 * i / len(cs), 360 * (i + 1) / len(cs), width=r * 0.22, fc=COR[c], ec="none", zorder=3))
    else:
        ax.add_patch(Wedge((x, y), r, 0, 360, width=r * 0.22, fc="#D0D0D0", ec="none", zorder=3))
    ax.text(x, y, row.rotulo, fontsize=6, ha="center", va="center", zorder=4)
mg = 0.08 * span
ax.set_xlim(xy[:, 0].min() - mg, xy[:, 0].max() + mg); ax.set_ylim(xy[:, 1].min() - mg, xy[:, 1].max() + mg)
ax.set_aspect("equal"); ax.axis("off")
cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.25, pad=0.01, location="left"); cb.set_label("log2FC (doente vs sadio)")
ax.legend(handles=[Patch(fc=COR[c], label=c) for c in CAT] + [Patch(fc="#D0D0D0", label="sem categoria")], loc="lower left", frameon=True, fontsize=9)
fig.savefig(f"{D}/subrede_{N}.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(f"{D}/subrede_{N}.pdf", bbox_inches="tight", facecolor="white")
print("Pronto.")
