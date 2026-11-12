# Protein-protein interaction network analysis (Methods 5.2.6): rede canais.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import numpy as np, pandas as pd, networkx as nx, gzip

D = "[diretório_transcriptoma]"
HITS = f"{D}/diamond_agave_canonical_degs_best_hits.tsv"
THR = 0.7; P = 0.041
OK = ["neighborhood", "fusion", "cooccurence", "coexpression", "experimental"]   # permitidos
ALL = OK + ["database", "textmining"]                                              # só para conferência

def comb(df, cols):
    prod = np.ones(len(df))
    for c in cols:
        prod *= 1 - np.clip((df[c].values - P) / (1 - P), 0, None)
    s = 1 - prod
    return s + P * (1 - s)

h = pd.read_csv(HITS, sep="\t", header=None, usecols=range(6),
                names=["transcrito", "ref", "evalue", "bitscore", "pident", "len"])
h = h.sort_values("bitscore", ascending=False).drop_duplicates("transcrito")
print(f"Transcritos canônicos com hit: {len(h)}")
rev = h.groupby("ref")["transcrito"].apply(list).to_dict()
refs = set(rev)

partes = []
for t in (4615, 906689):
    print(f"Lendo {t}...")
    for ch in pd.read_csv(f"{D}/string_db/{t}.protein.links.detailed.v12.0.txt.gz", sep=" ", chunksize=2_000_000):
        ch = ch[ch.protein1.isin(refs) & ch.protein2.isin(refs) & (ch.protein1 < ch.protein2)]
        if len(ch): partes.append(ch)
it = pd.concat(partes, ignore_index=True)
cols = ALL + ["combined_score"]
it[cols] = it[cols] / 1000
print(f"Pares de referência dentro do conjunto: {len(it)}")

it["s_all"] = comb(it, ALL)
print(f"Conferência (7 canais vs combined_score): diferença absoluta média = {(it.s_all - it.combined_score).abs().mean():.4f}")
it["s_ok"] = comb(it, OK)

for thr in (0.4, 0.7, 0.9):
    e = it[it.s_ok >= thr]
    print(f"score permitido >= {thr}: {len(e)} pares de referência, {len(set(e.protein1) | set(e.protein2))} proteínas de referência")

# expande para transcritos de Agave (hits podem ser compartilhados por vários transcritos)
e = it[it.s_ok >= THR]
linhas = []
for a, b, s in zip(e.protein1.values, e.protein2.values, e.s_ok.values):
    for x in rev[a]:
        for y in rev[b]:
            if x != y: linhas.append((x, y, a, b, round(s, 3)))
net = pd.DataFrame(linhas, columns=["agave_A", "agave_B", "ref_A", "ref_B", "score_permitido"]).drop_duplicates(["agave_A", "agave_B"])
G = nx.from_pandas_edgelist(net, "agave_A", "agave_B")
comps = sorted(nx.connected_components(G), key=len, reverse=True)
print(f"Rede em {THR}: {G.number_of_nodes()} nós, {G.number_of_edges()} arestas, {len(comps)} componentes; maior = {len(comps[0])} nós")
print(f"Componentes com >= 3 nós: {sum(len(c) >= 3 for c in comps)}")
net.to_csv(f"{D}/rede_canais_permitidos_thr{THR}.tsv", sep="\t", index=False)
print("Pronto.")
