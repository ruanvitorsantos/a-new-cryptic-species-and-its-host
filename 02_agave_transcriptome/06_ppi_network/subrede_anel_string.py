# Protein-protein interaction network analysis (Methods 5.2.6): subrede anel string.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import re, gzip, numpy as np, pandas as pd, networkx as nx
from functools import lru_cache
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Patch
from matplotlib.collections import LineCollection
from matplotlib.colors import TwoSlopeNorm
from matplotlib.cm import ScalarMappable

D = "[diretório_transcriptoma]"; N = 150; SEED = 42
nodes = pd.read_csv(f"{D}/subrede_{N}_nos.tsv", sep="\t", index_col=0)
edges = pd.read_csv(f"{D}/subrede_{N}_arestas.tsv", sep="\t")
full = pd.read_csv(f"{D}/rede_colapsada_nos.tsv", sep="\t", index_col=0)

par, cur = {}, None
for l in open(f"{D}/go-basic.obo"):
    l = l.strip()
    if l.startswith("["): cur = None
    if l == "[Term]": cur = "?"
    elif l.startswith("id: GO:") and cur == "?": cur = l[4:]; par[cur] = set()
    elif cur and cur != "?" and l.startswith("is_a: GO:"): par[cur].add(l[6:16])
    elif cur and cur != "?" and l.startswith("relationship: part_of GO:"): par[cur].add(l[22:32])
@lru_cache(None)
def anc(g): return frozenset({g}).union(*[anc(p) for p in par.get(g, ())]) if g in par else frozenset({g})

CAT = {"Resposta a estresse e defesa": ["GO:0006950", "GO:0006952"],
       "Controle de qualidade e degradação de proteínas": ["GO:0006457", "GO:0006508", "GO:0043161"],
       "Tradução e ribossomo": ["GO:0006412", "GO:0042254"],
       "Carboidratos e crescimento": ["GO:0005975", "GO:0040007"],
       "Replicação e reparo de DNA": ["GO:0006260", "GO:0006281"]}
COR = dict(zip(CAT, ["#D81B60", "#F2C230", "#8E5AB5", "#4CAF50", "#2B7BC9"]))

need = set(full.index); terms = {}
for t in (4615, 906689):
    with gzip.open(f"{D}/{t}.protein.enrichment.terms.v12.0.txt.gz", "rt") as f:
        for l in f:
            if "Biological Process" not in l: continue
            p = l.split("\t")
            if p[0] in need: terms.setdefault(p[0], set()).add(p[2])
def cats(r):
    A = set()
    for g in terms.get(r, ()): A |= anc(g)
    return [c for c, rs in CAT.items() if A & set(rs)]
full["cats"] = [cats(r) for r in full.index]
nodes["categorias"] = [";".join(cats(r)) for r in nodes.index]

for nome, x in (("Rede colapsada completa", full), ("Sub-rede", nodes)):
    dr = x.direction; cc = x["cats"] if "cats" in x else x.categorias.str.split(";")
    tab = {d: {c: 100 * np.mean([c in s for s in cc[dr == d]]) for c in CAT} for d in ("up", "down")}
    tab["n"] = {"up": int((dr == "up").sum()), "down": int((dr == "down").sum())}
    print(f"\n{nome} (% dos nós por categoria):"); print(pd.DataFrame(tab).round(1).to_string())
    print("Sem categoria:", int(sum(len(s) == 0 or s == [""] for s in cc)))

S = nx.from_pandas_edgelist(edges, "ref_A", "ref_B", edge_attr="score_permitido")
pos = nx.spring_layout(S, seed=SEED, k=1.6 / np.sqrt(len(S)), iterations=300, weight="score_permitido")
xy = np.array(list(pos.values())); span = max(np.ptp(xy[:, 0]), np.ptp(xy[:, 1])); r = 0.022 * span
top = set(nodes.sort_values("degree", ascending=False).head(25).index) | set(nodes.index[nodes.direction == "up"])
fig, ax = plt.subplots(figsize=(15, 15))
ax.add_collection(LineCollection([(pos[a], pos[b]) for a, b in S.edges()], colors="#8A8A8A",
                  linewidths=[0.4 + 2.2 * (S[a][b]["score_permitido"] - 0.7) / 0.3 for a, b in S.edges()], alpha=0.5, zorder=1))
norm = TwoSlopeNorm(vmin=-6, vcenter=0, vmax=6); cmap = plt.get_cmap("RdBu_r")
for n in S.nodes():
    x, y = pos[n]; row = nodes.loc[n]
    ax.add_patch(Circle((x, y), r * 0.8, fc=cmap(norm(np.clip(row.log2FoldChange, -6, 6))), ec="none", zorder=2))
    cs = [c for c in row.categorias.split(";") if c] or [None]
    for i, c in enumerate(cs):
        ax.add_patch(Wedge((x, y), r, 360 * i / len(cs), 360 * (i + 1) / len(cs), width=r * 0.22, fc=COR.get(c, "#D0D0D0"), ec="none", zorder=3))
    if n in top: ax.text(x, y, row.rotulo, fontsize=6, ha="center", va="center", zorder=4)
mg = 0.08 * span
ax.set_xlim(xy[:, 0].min() - mg, xy[:, 0].max() + mg); ax.set_ylim(xy[:, 1].min() - mg, xy[:, 1].max() + mg)
ax.set_aspect("equal"); ax.axis("off")
cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.25, pad=0.01, location="left"); cb.set_label("log2FC (doente vs sadio)")
ax.legend(handles=[Patch(fc=COR[c], label=c) for c in CAT] + [Patch(fc="#D0D0D0", label="sem categoria")], loc="lower left", fontsize=9)
fig.savefig(f"{D}/subrede_{N}_stringGO.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(f"{D}/subrede_{N}_stringGO.pdf", bbox_inches="tight", facecolor="white")
nodes.to_csv(f"{D}/subrede_{N}_nos_stringGO.tsv", sep="\t")
print("Pronto.")
