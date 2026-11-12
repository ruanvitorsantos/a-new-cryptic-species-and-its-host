# Protein-protein interaction network analysis (Methods 5.2.6): rede visual.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import re, gzip, numpy as np, pandas as pd, networkx as nx
from functools import lru_cache
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Patch
from matplotlib.collections import LineCollection
from matplotlib.colors import TwoSlopeNorm
from matplotlib.cm import ScalarMappable

D = "[diretório_transcriptoma]"
N = 150; SEED = 42; CAP_TRANS = 12; CAP_MOD = 30; FIGW = 30
full = pd.read_csv(f"{D}/rede_colapsada_nos.tsv", sep="\t", index_col=0)
edges = pd.read_csv(f"{D}/rede_colapsada_arestas.tsv", sep="\t")
G = nx.from_pandas_edgelist(edges, "ref_A", "ref_B", edge_attr="score_permitido")

# ---- GO (STRING) -> categorias
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
TR = "Tradução e ribossomo"
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

# ---- rótulos abreviados
info = pd.concat([pd.read_csv(f"{D}/string_db/{t}.protein.info.v12.0.txt.gz", sep="\t", usecols=[0, 3], names=["id", "ann"], header=0) for t in (4615, 906689)]).drop_duplicates("id").set_index("id")["ann"]
KEY = [("ribosomal", "RP"), ("heat shock", "HSP"), ("chaperon", "CHP"), ("ubiquitin", "UBQ"), ("proteasome", "PSM"),
       ("transcription factor", "TF"), ("elongation factor", "EF"), ("initiation factor", "IF"), ("kinase", "KIN"),
       ("phosphatase", "PP"), ("polymerase", "POL"), ("helicase", "HEL"), ("dehydrogenase", "DH"), ("transferase", "TRF"),
       ("synthase", "SYN"), ("hydrolase", "HYD"), ("oxidase", "OXI"), ("peroxidase", "POX"), ("cytochrome", "CYP"),
       ("v-type", "VATP"), ("histone", "HIS"), ("splicing", "SPL"), ("tubulin", "TUB"), ("actin", "ACT")]
def abrev(a):
    a = a.lower() if isinstance(a, str) else ""
    if not a or a.startswith("uncharacterized"): return "UNK"
    for k, v in KEY:
        if k in a: return v
    return re.split(r"[\s,;(]", a.strip())[0][:7].upper()
full["rotulo"] = [abrev(info.get(r, "")) for r in full.index]

# ---- seleção (critérios fixos, declarados antes)
full["escore"] = full[["degree", "betweenness"]].rank(pct=True).mean(axis=1)
sel, st, nmod = [], {"tr": 0}, {}
def add(r, force=False):
    if r in sel: return
    tr = TR in full.at[r, "cats"]; mo = full.at[r, "module"]
    if not force and ((tr and st["tr"] >= CAP_TRANS) or nmod.get(mo, 0) >= CAP_MOD): return
    sel.append(r); st["tr"] += int(tr); nmod[mo] = nmod.get(mo, 0) + 1
ups = list(full.index[full.direction == "up"])
for r in ups: add(r, force=True)
for r in ups:
    for v in sorted(G[r], key=lambda x: -full.at[x, "escore"])[:2]: add(v)
for r in full.sort_values("escore", ascending=False).index:
    if len(sel) >= N: break
    add(r)
S = G.subgraph(sel).copy()
for c in [c for c in nx.connected_components(S) if len(c) < 3]: S.remove_nodes_from(c)
m = full.loc[list(S.nodes())].copy()
comps = sorted(nx.connected_components(S), key=len, reverse=True)
print(f"Sub-rede: {S.number_of_nodes()} nós, {S.number_of_edges()} arestas, {len(comps)} componentes (maiores: {[len(c) for c in comps[:8]]})")
print(f"Up: {(m.direction == 'up').sum()} | Down: {(m.direction == 'down').sum()} | nós de tradução: {sum(TR in c for c in m.cats)}")

# ---- layout por componente + empacotamento
def separar(P, dmin=4.6, it=800):
    P = P.copy(); I = np.eye(len(P), dtype=bool)
    for _ in range(it):
        d = P[:, None, :] - P[None, :, :]; dist = np.linalg.norm(d, axis=2) + 1e-9
        close = (dist < dmin) & ~I
        if not close.any(): break
        P += ((d / dist[..., None]) * ((dmin - dist) / 2)[..., None] * close[..., None]).sum(axis=1)
    return P
blocks = []
for c in comps:
    H = S.subgraph(c); n = len(H); nodes = list(H.nodes())
    p = nx.spring_layout(H, k=3.0 / np.sqrt(n), iterations=600, seed=SEED, scale=1.0, weight="score_permitido")
    P = separar(np.array([p[x] for x in nodes]) * (4.2 * np.sqrt(n) + 2.0))
    P -= P.min(axis=0); blocks.append((nodes, P, P[:, 0].max() + 3, P[:, 1].max() + 3))
Wt = max(np.sqrt(sum(w * h for _, _, w, h in blocks)) * 1.3, max(b[2] for b in blocks))
pos, x, y, rowh = {}, 0, 0, 0
for nodes, P, w, h in blocks:
    if x + w > Wt: x, y, rowh = 0, y + rowh, 0
    for n_, q in zip(nodes, P): pos[n_] = (q[0] + x + 1.5, -(q[1] + y + 1.5))
    x += w; rowh = max(rowh, h)
Ht = y + rowh

# ---- figura
xy = np.array(list(pos.values())); mg = 4
W = xy[:, 0].max() - xy[:, 0].min() + 2 * mg; Hh = xy[:, 1].max() - xy[:, 1].min() + 2 * mg
fig, ax = plt.subplots(figsize=(FIGW, max(FIGW * Hh / W, 6)))
ppu = FIGW * 72 / W; fs = float(np.clip(0.5 * ppu, 4, 12)); r = 1.0
ax.add_collection(LineCollection([(pos[a], pos[b]) for a, b in S.edges()], colors="#8A8A8A", alpha=0.55, zorder=1,
                  linewidths=[0.5 + 2.0 * (S[a][b]["score_permitido"] - 0.7) / 0.3 for a, b in S.edges()]))
_tp = lambda a, b: "-".join(sorted([m.at[a, "direction"], m.at[b, "direction"]]))
for _t, _c, _w in (("up-up", "#D7191C", 2.6), ("down-up", "#7B3294", 2.2)):
    _sg = [(pos[a], pos[b]) for a, b in S.edges() if _tp(a, b) == _t]
    if _sg: ax.add_collection(LineCollection(_sg, colors=_c, linewidths=_w, alpha=0.9, zorder=1.5))
norm = TwoSlopeNorm(vmin=-6, vcenter=0, vmax=6); cmap = plt.get_cmap("RdBu_r")
for n_ in S.nodes():
    cx, cy = pos[n_]; row = m.loc[n_]
    if row.direction == "up": ax.add_patch(Circle((cx, cy), r * 1.4, fc="none", ec="#E31A1C", lw=2.4, zorder=6))
    ax.add_patch(Circle((cx, cy), r * 0.78, fc=cmap(norm(np.clip(row.log2FoldChange, -6, 6))), ec="none", zorder=2))
    cs = row.cats
    if not cs: ax.add_patch(Circle((cx, cy), r * 0.9, fc="none", ec="#B0B0B0", lw=0.8, zorder=3))
    for i, c in enumerate(cs):
        ax.add_patch(Wedge((cx, cy), r, 360 * i / len(cs), 360 * (i + 1) / len(cs), width=r * 0.24, fc=COR.get(c, "#D0D0D0"), ec="none", zorder=3))
    lg = len(row.rotulo) > 5
    ax.text(cx, cy - (r * 1.25 if lg else 0), row.rotulo, fontsize=fs, ha="center", va="top" if lg else "center", zorder=4,
            bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.4) if lg else None)
ax.set_xlim(xy[:, 0].min() - mg, xy[:, 0].max() + mg); ax.set_ylim(xy[:, 1].min() - mg, xy[:, 1].max() + mg)
ax.set_aspect("equal"); ax.axis("off")
cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.25, pad=0.01, location="left"); cb.set_label("log2FC (doente vs sadio)")
ax.legend(handles=[Patch(fc=COR[c], label=c) for c in CAT] + [Patch(fc="#D0D0D0", label="outras / sem categoria")], loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=10)
from matplotlib.lines import Line2D
from collections import Counter
print("Tipos de aresta na figura:", dict(Counter("-".join(sorted([m.at[a, "direction"], m.at[b, "direction"]])) for a, b in S.edges())))
_l1 = ax.get_legend(); ax.add_artist(_l1)
ax.legend(handles=[Line2D([0], [0], marker="o", ls="", mfc="white", mec="#E31A1C", mew=2.4, ms=11, label="gene UP (halo vermelho)"),
                   Line2D([0], [0], marker="o", ls="", mfc="white", mec="#B0B0B0", ms=11, label="gene DOWN (sem halo)"),
                   Line2D([0], [0], color="#D7191C", lw=2.6, label="interação UP-UP"),
                   Line2D([0], [0], color="#7B3294", lw=2.2, label="interação UP-DOWN"),
                   Line2D([0], [0], color="#8A8A8A", lw=1.2, label="interação DOWN-DOWN")],
          loc="upper left", bbox_to_anchor=(1.01, 0.95), fontsize=10)
fig.savefig(f"{D}/rede_visual_updown.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(f"{D}/rede_visual_updown.pdf", bbox_inches="tight", facecolor="white")
m["x"] = [pos[i][0] for i in m.index]; m["y"] = [pos[i][1] for i in m.index]
m.assign(cats=m.cats.apply(";".join)).to_csv(f"{D}/rede_visual_nos.tsv", sep="\t")
nx.to_pandas_edgelist(S).rename(columns={"source": "ref_A", "target": "ref_B"}).to_csv(f"{D}/rede_visual_arestas.tsv", sep="\t", index=False)
print("Pronto.")
