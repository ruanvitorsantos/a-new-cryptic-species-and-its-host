# Protein-protein interaction network analysis (Methods 5.2.6): comparar redes.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import os, numpy as np, pandas as pd, networkx as nx
D = "[diretório_transcriptoma]"
DEG = "[diretório_deseq2_results]/DEGs_todos_shrink.tsv"
OUT = f"{D}/comparacao5"; os.makedirs(OUT, exist_ok=True)
THR = 0.7; P = 0.041; SEED = 42; MINPID = float(os.environ.get("MINPID", "0"))
OK = ["neighborhood", "fusion", "cooccurence", "coexpression", "experimental"]; ALL = OK + ["database", "textmining"]
tid = lambda s: s.str.replace(r"\.p\d+$", "", regex=True)
def comb(df, cols):
    prod = np.ones(len(df))
    for c in cols: prod *= 1 - np.clip((df[c].values - P) / (1 - P), 0, None)
    s = 1 - prod; return s + P * (1 - s)
dg = pd.read_csv(DEG, sep="\t").set_index("gene")

h5 = pd.read_csv(f"{D}/ref5/hits5_all.tsv", sep="\t", header=None, names=["q", "s", "e", "b", "p", "l"])
h5["tax"] = h5.s.str.split(".").str[0]; h5 = h5[h5.p >= MINPID]
ad = pd.read_csv(f"{D}/diamond_agave_canonical_degs_best_hits.tsv", sep="\t", header=None, usecols=[0, 1, 4], names=["q", "s", "p"])
maps = {"ananas+dendrobium": ad[ad.p >= MINPID]}
for t, n in [("3702", "arabidopsis"), ("4577", "milho"), ("39947", "arroz"), ("4558", "sorgo"), ("214687", "banana")]:
    maps[n] = h5[h5.tax == t].sort_values("b", ascending=False).drop_duplicates("q")[["q", "s", "p"]]
maps["global5"] = h5.sort_values("b", ascending=False).drop_duplicates("q")[["q", "s", "p"]]
refs_tax = {}
for m in maps.values():
    for s in m.s.unique(): refs_tax.setdefault(s.split(".")[0], set()).add(s)

edges = {}
for t, refs in refs_tax.items():
    d = "string_db" if t in ("4615", "906689") else "string_db5"
    parts, dif, cnt = [], 0.0, 0
    for ch in pd.read_csv(f"{D}/{d}/{t}.protein.links.detailed.v12.0.txt.gz", sep=" ", chunksize=2_000_000):
        ch = ch[ch.protein1.isin(refs) & ch.protein2.isin(refs) & (ch.protein1 < ch.protein2)].copy()
        if ch.empty: continue
        ch[ALL + ["combined_score"]] = ch[ALL + ["combined_score"]] / 1000
        dif += (comb(ch, ALL) - ch.combined_score).abs().sum(); cnt += len(ch)
        ch["s"] = comb(ch, OK); parts.append(ch.loc[ch.s >= THR, ["protein1", "protein2", "s"]])
    edges[t] = pd.concat(parts) if parts else pd.DataFrame(columns=["protein1", "protein2", "s"])
    print(f"{t}: pares verificados = {cnt}, diferença média do recálculo = {dif / max(cnt, 1):.4f}, arestas >= {THR}: {len(edges[t])}")
E_all = pd.concat(edges.values())

res, nodes_t, edges_t, mod_t = [], {}, {}, {}
for nome, m in maps.items():
    m = m.copy(); m["t"] = tid(m.q); m = m.drop_duplicates("t"); refs = set(m.s)
    E = E_all[E_all.protein1.isin(refs) & E_all.protein2.isin(refs)]
    G = nx.Graph(); G.add_weighted_edges_from(zip(E.protein1, E.protein2, E.s))
    for c in [c for c in nx.connected_components(G) if len(c) < 3]: G.remove_nodes_from(c)
    if G.number_of_nodes() == 0: print(nome, "vazia"); continue
    com = nx.community.louvain_communities(G, seed=SEED); modref = {n: i + 1 for i, c in enumerate(sorted(com, key=len, reverse=True)) for n in c}
    tr = m[m.s.isin(G.nodes())].copy()
    tr["lfc"] = dg["log2FoldChange"].reindex(tr.t).values; tr["padj"] = dg["padj"].reindex(tr.t).values
    tr["modulo"] = tr.s.map(modref)
    rep = tr.sort_values("padj").drop_duplicates("s")
    byref = tr.groupby("s").t.apply(list)
    et = set()
    for a, b in G.edges():
        for x in byref[a]:
            for y in byref[b]:
                if x != y: et.add((x, y) if x < y else (y, x))
    nodes_t[nome] = set(tr.t); edges_t[nome] = et; mod_t[nome] = dict(zip(tr.t, tr.modulo))
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    res.append(dict(rede=nome, transcritos_com_hit=len(m), pident_mediana=round(m.p.median(), 1), refs_distintas=m.s.nunique(),
                    nos_ref=G.number_of_nodes(), arestas_ref=G.number_of_edges(), componentes=len(comps), maior_comp=len(comps[0]),
                    transcritos_na_rede=len(tr), up_transc=int((tr.lfc > 0).sum()), down_transc=int((tr.lfc < 0).sum()),
                    up_ref=int((rep.lfc > 0).sum()), down_ref=int((rep.lfc < 0).sum()),
                    densidade=round(nx.density(G), 4), modularidade=round(nx.community.modularity(G, com), 3),
                    modulos_10=sum(len(c) >= 10 for c in com)))
    deg = pd.Series(dict(G.degree()))
    rep.assign(degree=rep.s.map(deg), direcao=np.where(rep.lfc > 0, "up", "down")).to_csv(f"{OUT}/rede_{nome}_nos.tsv", sep="\t", index=False)
    pd.DataFrame([(a, b, d["weight"]) for a, b, d in G.edges(data=True)], columns=["ref_A", "ref_B", "score_permitido"]).to_csv(f"{OUT}/rede_{nome}_arestas.tsv", sep="\t", index=False)

R = pd.DataFrame(res).set_index("rede"); R.to_csv(f"{OUT}/resumo_redes.tsv", sep="\t")
print("\n" + R.T.to_string())
nm = list(R.index)
jac = lambda a, b: len(a & b) / len(a | b) if (a | b) else np.nan
def ari(a, b):
    ct = pd.crosstab(a, b).values.astype(float); n = ct.sum(); c2 = lambda x: x * (x - 1) / 2
    s, sa, sb = c2(ct).sum(), c2(ct.sum(1)).sum(), c2(ct.sum(0)).sum(); ex = sa * sb / c2(n); mx = (sa + sb) / 2
    return (s - ex) / (mx - ex) if mx != ex else np.nan
JN = pd.DataFrame([[jac(nodes_t[a], nodes_t[b]) for b in nm] for a in nm], index=nm, columns=nm)
JE = pd.DataFrame([[jac(edges_t[a], edges_t[b]) for b in nm] for a in nm], index=nm, columns=nm)
def ari_par(a, b):
    sh = sorted(nodes_t[a] & nodes_t[b])
    return ari(pd.Series([mod_t[a][x] for x in sh]), pd.Series([mod_t[b][x] for x in sh])) if len(sh) >= 10 else np.nan
AR = pd.DataFrame([[ari_par(a, b) for b in nm] for a in nm], index=nm, columns=nm)
for t, nmx in ((JN, "jaccard_nos"), (JE, "jaccard_arestas"), (AR, "ari_modulos")):
    t.round(3).to_csv(f"{OUT}/{nmx}.tsv", sep="\t"); print(f"\n{nmx}:"); print(t.round(2).to_string())
print("\nPronto:", OUT)
