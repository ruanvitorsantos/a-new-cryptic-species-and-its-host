# Protein-protein interaction network analysis (Methods 5.2.6): rede abacaxi.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import os, re, gzip, numpy as np, pandas as pd, networkx as nx
from functools import lru_cache
from collections import Counter
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Patch
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from matplotlib.colors import TwoSlopeNorm
from matplotlib.cm import ScalarMappable

D = "[diretório_transcriptoma]"
DEG = "[diretório_deseq2_results]/DEGs_todos_shrink.tsv"
OUT = f"{D}/rede_abacaxi"; os.makedirs(OUT, exist_ok=True)
THR = 0.7; P = 0.041; NSEL = 150; SEED = 42; CAP_TRANS = 12; CAP_MOD = 30; FIGW = 30
OK = ["neighborhood", "fusion", "cooccurence", "coexpression", "experimental"]
tid = lambda s: re.sub(r"\.p\d+$", "", s)
def comb(df, cols):
    prod = np.ones(len(df))
    for c in cols: prod *= 1 - np.clip((df[c].values - P) / (1 - P), 0, None)
    s = 1 - prod; return s + P * (1 - s)
dg = pd.read_csv(DEG, sep="\t").set_index("gene")

# ---- transcritos com hit no abacaxi (canônicos + recuperados)
h = pd.read_csv(f"{D}/ref7/hits7_all_plus.tsv", sep="\t", header=None, names=["q", "s", "e", "b", "p", "l"])
h["tax"] = h.s.str.split(".").str[0]; h["t"] = h.q.map(tid)
rec = set(pd.read_csv(f"{D}/ref7/recuperados_info.tsv", sep="\t", index_col=0).index)
M = h[h.tax == "4615"].sort_values("b", ascending=False).drop_duplicates("t")[["t", "s", "b", "p"]].copy()
M["lfc"] = dg.log2FoldChange.reindex(M.t).values; M["padj"] = dg.padj.reindex(M.t).values; M["rec"] = M.t.isin(rec)
ara = h[h.tax == "3702"].sort_values("b", ascending=False).drop_duplicates("t").set_index("t").s
print(f"Transcritos com hit no abacaxi: {len(M)} (up {(M.lfc > 0).sum()}, down {(M.lfc < 0).sum()}) | recuperados: {int(M.rec.sum())} (up {int((M.rec & (M.lfc > 0)).sum())})")

# ---- arestas do abacaxi, só canais permitidos
refs = set(M.s); parts = []
for ch in pd.read_csv(f"{D}/string_db/4615.protein.links.detailed.v12.0.txt.gz", sep=" ", chunksize=2_000_000):
    ch = ch[ch.protein1.isin(refs) & ch.protein2.isin(refs) & (ch.protein1 < ch.protein2)].copy()
    if ch.empty: continue
    ch[OK] = ch[OK] / 1000; ch["score"] = comb(ch, OK)
    parts.append(ch.loc[ch.score >= THR, ["protein1", "protein2", "score"]])
E = pd.concat(parts)
G = nx.Graph(); G.add_weighted_edges_from(zip(E.protein1, E.protein2, E.score), weight="score")
for c in [c for c in nx.connected_components(G) if len(c) < 3]: G.remove_nodes_from(c)
R = M.sort_values("padj").drop_duplicates("s").set_index("s").loc[list(G.nodes())].copy()
R["n_transcritos"] = M.groupby("s").size().reindex(R.index)
R["degree"] = pd.Series(dict(G.degree())); R["betweenness"] = pd.Series(nx.betweenness_centrality(G))
com = sorted(nx.community.louvain_communities(G, seed=SEED), key=len, reverse=True)
R["module"] = pd.Series({n: i + 1 for i, c in enumerate(com) for n in c})
R["log2FoldChange"] = R.lfc; R["direction"] = np.where(R.lfc > 0, "up", "down")
print(f"Rede completa do abacaxi: {G.number_of_nodes()} nós, {G.number_of_edges()} arestas, {len(com)} módulos | up {(R.direction == 'up').sum()}, down {(R.direction == 'down').sum()} | recuperados: up {int((R.rec & (R.direction == 'up')).sum())}, down {int((R.rec & (R.direction == 'down')).sum())}")

# ---- anotação: rótulo, organela, categorias GO
info = pd.read_csv(f"{D}/string_db/4615.protein.info.v12.0.txt.gz", sep="\t", usecols=[0, 3], names=["id", "ann"], header=0).set_index("id")["ann"]
isym = pd.read_csv(f"{D}/string_db5/3702.protein.info.v12.0.txt.gz", sep="\t", usecols=[0, 1], names=["id", "sym"], header=0).set_index("id")["sym"]
R["organela"] = [bool(re.match(r"^(atp|nad|ndh|rps|rpl|psa|psb|pet|rbc|cox|cob|ccm|mat|rpo|ycf)", str(isym.get(ara.get(t, ""), "")))) for t in R.t]
KEY = [("ribosomal", "RP"), ("heat shock", "HSP"), ("chaperon", "CHP"), ("ubiquitin", "UBQ"), ("proteasome", "PSM"), ("transcription factor", "TF"),
       ("elongation factor", "EF"), ("initiation factor", "IF"), ("kinase", "KIN"), ("phosphatase", "PP"), ("polymerase", "POL"), ("helicase", "HEL"),
       ("dehydrogenase", "DH"), ("transferase", "TRF"), ("synthase", "SYN"), ("hydrolase", "HYD"), ("oxidase", "OXI"), ("peroxidase", "POX"),
       ("cytochrome", "CYP"), ("v-type", "VATP"), ("histone", "HIS"), ("splicing", "SPL"), ("tubulin", "TUB"), ("actin", "ACT")]
GEN = {"PROTEIN", "PUTATIV", "PUTATIVE", "SMALL", "MATRIN-", "PROBABLE", "UNCHARA", "UNCHARACTERIZED"}
def abrev(a):
    a = a.lower() if isinstance(a, str) else ""
    if not a or a.startswith("uncharacterized"): return "UNK"
    for k, v in KEY:
        if k in a: return v
    w = re.split(r"[\s,;(]", a.strip())[0][:7].upper()
    return "UNK" if w in GEN else w
R["rotulo"] = [abrev(info.get(r, "")) for r in R.index]
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
CAT = {"Resposta a estresse e defesa": ["GO:0006950", "GO:0006952"], "Controle de qualidade e degradação de proteínas": ["GO:0006457", "GO:0006508", "GO:0043161"],
       "Tradução e ribossomo": ["GO:0006412", "GO:0042254"], "Carboidratos e crescimento": ["GO:0005975", "GO:0040007"], "Replicação e reparo de DNA": ["GO:0006260", "GO:0006281"]}
COR = dict(zip(CAT, ["#D81B60", "#F2C230", "#8E5AB5", "#4CAF50", "#2B7BC9"])); TR = "Tradução e ribossomo"
need = set(R.index); terms = {}
with gzip.open(f"{D}/4615.protein.enrichment.terms.v12.0.txt.gz", "rt") as f:
    for l in f:
        if "Biological Process" not in l: continue
        p = l.split("\t")
        if p[0] in need: terms.setdefault(p[0], set()).add(p[2])
def cats(r):
    A = set()
    for g in terms.get(r, ()): A |= anc(g)
    return [c for c, rs in CAT.items() if A & set(rs)]
R["cats"] = [cats(r) for r in R.index]

# ---- seleção da sub-rede (mesma regra das figuras anteriores)
R["escore"] = R[["degree", "betweenness"]].rank(pct=True).mean(axis=1)
sel, st, nmod = [], {"tr": 0}, {}
def add(r, force=False):
    if r in sel: return
    tr = TR in R.at[r, "cats"]; mo = R.at[r, "module"]
    if not force and ((tr and st["tr"] >= CAP_TRANS) or nmod.get(mo, 0) >= CAP_MOD): return
    sel.append(r); st["tr"] += int(tr); nmod[mo] = nmod.get(mo, 0) + 1
ups = list(R.index[R.direction == "up"])
for r in ups: add(r, force=True)
for r in ups:
    for v in sorted(G[r], key=lambda x: -R.at[x, "escore"])[:2]: add(v)
for r in R.sort_values("escore", ascending=False).index:
    if len(sel) >= NSEL: break
    add(r)
S = G.subgraph(sel).copy()
for c in [c for c in nx.connected_components(S) if len(c) < 3]: S.remove_nodes_from(c)
m = R.loc[list(S.nodes())].copy()
comps = sorted(nx.connected_components(S), key=len, reverse=True)
tp = Counter("-".join(sorted([m.at[a, "direction"], m.at[b, "direction"]])) for a, b in S.edges())
print(f"Figura: {S.number_of_nodes()} nós, {S.number_of_edges()} arestas, {len(comps)} componentes (maiores: {[len(c) for c in comps[:8]]})")
print(f"Up {(m.direction == 'up').sum()} de {len(ups)} | Down {(m.direction == 'down').sum()} | organela {int(m.organela.sum())} (up {int((m.organela & (m.direction == 'up')).sum())}) | recuperados {int(m.rec.sum())} | tradução {sum(TR in c for c in m.cats)}")
print("Arestas por tipo:", dict(tp))

# ---- layout por componente
def separar(Pp, dmin=4.6, it=800):
    Pp = Pp.copy(); I = np.eye(len(Pp), dtype=bool)
    for _ in range(it):
        d = Pp[:, None, :] - Pp[None, :, :]; dist = np.linalg.norm(d, axis=2) + 1e-9
        close = (dist < dmin) & ~I
        if not close.any(): break
        Pp += ((d / dist[..., None]) * ((dmin - dist) / 2)[..., None] * close[..., None]).sum(axis=1)
    return Pp
blocks = []
for c in comps:
    H = S.subgraph(c); n = len(H); nodes = list(H.nodes())
    p = nx.spring_layout(H, k=3.0 / np.sqrt(n), iterations=600, seed=SEED, scale=1.0, weight="score")
    Pp = separar(np.array([p[x] for x in nodes]) * (4.2 * np.sqrt(n) + 2.0))
    Pp -= Pp.min(axis=0); blocks.append((nodes, Pp, Pp[:, 0].max() + 3, Pp[:, 1].max() + 3))
Wt = max(np.sqrt(sum(w * hh for _, _, w, hh in blocks)) * 1.3, max(b[2] for b in blocks))
pos, x, y, rowh = {}, 0, 0, 0
for nodes, Pp, w, hh in blocks:
    if x + w > Wt: x, y, rowh = 0, y + rowh, 0
    for n_, q in zip(nodes, Pp): pos[n_] = (q[0] + x + 1.5, -(q[1] + y + 1.5))
    x += w; rowh = max(rowh, hh)

# ---- figura
xy = np.array(list(pos.values())); mg = 4
W_ = xy[:, 0].max() - xy[:, 0].min() + 2 * mg; Hh = xy[:, 1].max() - xy[:, 1].min() + 2 * mg
fig, ax = plt.subplots(figsize=(FIGW, max(FIGW * Hh / W_, 6)))
fs = float(np.clip(0.5 * FIGW * 72 / W_, 4, 12)); r = 1.0
ax.add_collection(LineCollection([(pos[a], pos[b]) for a, b in S.edges()], colors="#8A8A8A", alpha=0.55, zorder=1,
                  linewidths=[0.5 + 2.0 * (S[a][b]["score"] - 0.7) / 0.3 for a, b in S.edges()]))
for t_, c_, w_ in (("up-up", "#D7191C", 2.6), ("down-up", "#7B3294", 2.2)):
    sg = [(pos[a], pos[b]) for a, b in S.edges() if "-".join(sorted([m.at[a, "direction"], m.at[b, "direction"]])) == t_]
    if sg: ax.add_collection(LineCollection(sg, colors=c_, linewidths=w_, alpha=0.9, zorder=1.5))
norm = TwoSlopeNorm(vmin=-6, vcenter=0, vmax=6); cmap = plt.get_cmap("RdBu_r")
for n_ in S.nodes():
    cx, cy = pos[n_]; row = m.loc[n_]
    if row.direction == "up": ax.add_patch(Circle((cx, cy), r * 1.4, fc="none", ec="#E31A1C", lw=2.4, zorder=6))
    if row.organela: ax.add_patch(Circle((cx, cy), r * 1.12, fc="none", ec="black", lw=1.6, ls="--", zorder=5))
    ax.add_patch(Circle((cx, cy), r * 0.78, fc=cmap(norm(np.clip(row.log2FoldChange, -6, 6))), ec="#B0B0B0", lw=0.6, zorder=2))
    cs = row.cats
    for i, c in enumerate(cs):
        ax.add_patch(Wedge((cx, cy), r, 360 * i / len(cs), 360 * (i + 1) / len(cs), width=r * 0.24, fc=COR[c], ec="none", zorder=3))
    if row.rec: ax.plot(cx + 0.95 * r, cy + 0.95 * r, marker="s", ms=5, color="black", zorder=7)
    lg = len(row.rotulo) > 5
    ax.text(cx, cy - (r * 1.25 if lg else 0), row.rotulo, fontsize=fs, ha="center", va="top" if lg else "center", zorder=4,
            bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.4) if lg else None)
ax.set_xlim(xy[:, 0].min() - mg, xy[:, 0].max() + mg); ax.set_ylim(xy[:, 1].min() - mg, xy[:, 1].max() + mg)
ax.set_aspect("equal"); ax.axis("off")
cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.25, pad=0.01, location="left"); cb.set_label("log2FC (doente vs sadio)")
leg1 = ax.legend(handles=[Patch(fc=COR[c], label=c) for c in CAT] + [Patch(fc="white", ec="#B0B0B0", label="outras / sem categoria"),
    Line2D([0], [0], marker="o", ls="", mfc="white", mec="black", ms=11, label="gene de organela (contorno tracejado)"),
    Line2D([0], [0], marker="s", ls="", color="black", ms=6, label="recuperado por blastx (sem ORF predito)")],
    loc="center left", bbox_to_anchor=(1.01, 0.66), fontsize=10)
ax.add_artist(leg1)
ax.legend(handles=[Line2D([0], [0], marker="o", ls="", mfc="white", mec="#E31A1C", mew=2.4, ms=11, label="gene UP (halo vermelho)"),
    Line2D([0], [0], marker="o", ls="", mfc="white", mec="#B0B0B0", ms=11, label="gene DOWN (sem halo)"),
    Line2D([0], [0], color="#D7191C", lw=2.6, label="interação UP-UP"), Line2D([0], [0], color="#7B3294", lw=2.2, label="interação UP-DOWN"),
    Line2D([0], [0], color="#8A8A8A", lw=1.2, label="interação DOWN-DOWN")], loc="center left", bbox_to_anchor=(1.01, 0.28), fontsize=10)
for ext in ("png", "pdf"): fig.savefig(f"{OUT}/rede_abacaxi.{ext}", dpi=300, bbox_inches="tight", facecolor="white")

# ---- arquivos para o Cytoscape
m["x"] = [pos[i][0] for i in m.index]; m["y"] = [pos[i][1] for i in m.index]
cc = {"c_estresse": "Resposta a estresse", "c_proteinas": "Controle de qualidade", "c_traducao": "Tradução", "c_carboidratos": "Carboidratos", "c_dna": "Replicação e reparo"}
cs_ = m.cats.apply(";".join)
for k, v in cc.items(): m[k] = cs_.str.contains(v, regex=False).astype(int)
m["c_outras"] = (m[list(cc)].sum(axis=1) == 0).astype(int); m.index.name = "id"
m[["t", "rotulo", "direction", "log2FoldChange", "padj", "module", "organela", "rec", "n_transcritos", "x", "y"] + list(cc) + ["c_outras"]].to_csv(f"{OUT}/nos.tsv", sep="\t")
ed = pd.DataFrame([(a, b, d["score"], "-".join(sorted([m.at[a, "direction"], m.at[b, "direction"]]))) for a, b, d in S.edges(data=True)], columns=["source", "target", "peso", "tipo"])
ed.to_csv(f"{OUT}/arestas.tsv", sep="\t", index=False)
print("Pronto:", OUT)
Rn = R.copy(); Rn.index.name = "id"
cs_f = Rn.cats.apply(";".join)
for k, v in cc.items(): Rn[k] = cs_f.str.contains(v, regex=False).astype(int)
Rn["c_outras"] = (Rn[list(cc)].sum(axis=1) == 0).astype(int)
Rn["na_figura"] = Rn.index.isin(m.index)
Rn[["t", "rotulo", "direction", "log2FoldChange", "padj", "module", "degree", "betweenness", "organela", "rec", "n_transcritos", "na_figura"] + list(cc) + ["c_outras"]].to_csv(f"{OUT}/nos_completa.tsv", sep="\t")
ef = pd.DataFrame([(a, b, d["score"], "-".join(sorted([R.at[a, "direction"], R.at[b, "direction"]]))) for a, b, d in G.edges(data=True)], columns=["source", "target", "peso", "tipo"])
ef.to_csv(f"{OUT}/arestas_completa.tsv", sep="\t", index=False)
print("Rede completa exportada:", len(Rn), "nós,", len(ef), "arestas")
