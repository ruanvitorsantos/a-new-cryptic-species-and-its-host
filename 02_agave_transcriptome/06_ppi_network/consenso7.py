# Protein-protein interaction network analysis (Methods 5.2.6): consenso7.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import os, numpy as np, pandas as pd, networkx as nx
from collections import defaultdict
D = "[diretório_transcriptoma]"
DEG = "[diretório_deseq2_results]/DEGs_todos_shrink.tsv"
OUT = f"{D}/consenso7"; os.makedirs(OUT, exist_ok=True)
THR = 0.7; P = 0.041; KSEL = int(os.environ.get("KSEL", "3")); SEED = 42
OK = ["neighborhood", "fusion", "cooccurence", "coexpression", "experimental"]
TAX = ["4615", "906689", "3702", "4577", "39947", "4558", "214687"]
NOME = dict(zip(TAX, ["abacaxi", "dendrobium", "arabidopsis", "milho", "arroz", "sorgo", "banana"]))
tid = lambda s: s.str.replace(r"\.p\d+$", "", regex=True)
def comb(df, cols):
    prod = np.ones(len(df))
    for c in cols: prod *= 1 - np.clip((df[c].values - P) / (1 - P), 0, None)
    s = 1 - prod; return s + P * (1 - s)
dg = pd.read_csv(DEG, sep="\t").set_index("gene")

h = pd.read_csv(f"{D}/ref7/hits7_all.tsv", sep="\t", header=None, names=["q", "s", "e", "b", "p", "l"])
h["tax"] = h.s.str.split(".").str[0]
best = h.sort_values("b", ascending=False).drop_duplicates(["q", "tax"])
W = best.pivot(index="q", columns="tax", values="s").reindex(columns=TAX)
W["bmax"] = best.groupby("q").b.max(); W["t"] = tid(pd.Series(W.index, index=W.index))
W = W.sort_values("bmax", ascending=False).drop_duplicates("t"); W.index = W.t
W["lfc"] = dg["log2FoldChange"].reindex(W.index).values; W["padj"] = dg["padj"].reindex(W.index).values
W["chave"] = W[TAX].fillna("-").agg("|".join, axis=1)
R = W.sort_values("padj").drop_duplicates("chave")
print(f"Transcritos com hit: {len(W)} | nós após colapsar idênticos nas 7 referências: {len(R)}")
HIT = {t: R.loc[t, TAX].notna().values for t in R.index}

sup = defaultdict(set)
for t in TAX:
    refs = set(R[t].dropna()); byref = R.groupby(t).apply(lambda g: list(g.index)).to_dict()
    n = 0
    for ch in pd.read_csv(f"{D}/{'string_db' if t in ('4615', '906689') else 'string_db5'}/{t}.protein.links.detailed.v12.0.txt.gz", sep=" ", chunksize=2_000_000):
        ch = ch[ch.protein1.isin(refs) & ch.protein2.isin(refs) & (ch.protein1 < ch.protein2)].copy()
        if ch.empty: continue
        ch[OK] = ch[OK] / 1000; ch = ch[comb(ch, OK) >= THR]
        for a, b in zip(ch.protein1.values, ch.protein2.values):
            for x in byref[a]:
                for y in byref[b]:
                    if x != y: sup[(x, y) if x < y else (y, x)].add(t); n += 1
    print(f"{NOME[t]}: arestas de referência >= {THR}: {n}")
E = pd.DataFrame([(a, b, len(s), ";".join(NOME[x] for x in TAX if x in s), int((HIT[a] & HIT[b]).sum())) for (a, b), s in sup.items()],
                 columns=["A", "B", "K", "referencias", "possivel"])
E.to_csv(f"{OUT}/arestas_consenso_todas.tsv", sep="\t", index=False)

def rede(k):
    G = nx.from_pandas_edgelist(E[E.K >= k], "A", "B", edge_attr="K")
    for c in [c for c in nx.connected_components(G) if len(c) < 3]: G.remove_nodes_from(c)
    return G
lin = []
for k in range(1, 8):
    G = rede(k)
    if G.number_of_nodes() == 0: lin.append(dict(K=k, arestas=0)); continue
    cs = sorted(nx.connected_components(G), key=len, reverse=True); l = R.loc[list(G.nodes()), "lfc"]
    lin.append(dict(K=k, nos=G.number_of_nodes(), arestas=G.number_of_edges(), componentes=len(cs), maior=len(cs[0]),
                    up=int((l > 0).sum()), down=int((l < 0).sum()), pct_up=round(100 * (l > 0).mean(), 1), densidade=round(nx.density(G), 4)))
print("\nRede de consenso por K (K = nº de referências que sustentam a aresta):"); print(pd.DataFrame(lin).to_string(index=False))

G = rede(KSEL)
if G.number_of_nodes():
    com = sorted(nx.community.louvain_communities(G, seed=SEED), key=len, reverse=True)
    mod = {n: i + 1 for i, c in enumerate(com) for n in c}
    N = R.loc[list(G.nodes())].copy()
    N["degree"] = pd.Series(dict(G.degree())); N["betweenness"] = pd.Series(nx.betweenness_centrality(G))
    N["module"] = pd.Series(mod); N["direction"] = np.where(N.lfc > 0, "up", "down")
    N.rename(columns=NOME).drop(columns=["chave", "t", "bmax"]).to_csv(f"{OUT}/nos_consenso_K{KSEL}.tsv", sep="\t")
    nx.to_pandas_edgelist(G).rename(columns={"source": "A", "target": "B"}).to_csv(f"{OUT}/arestas_consenso_K{KSEL}.tsv", sep="\t", index=False)
    print(f"\nRede salva: K >= {KSEL}, {G.number_of_nodes()} nós, {G.number_of_edges()} arestas")
print("Pronto:", OUT)
