# Protein-protein interaction network analysis (Methods 5.2.6): colapsar rede.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import numpy as np, pandas as pd, networkx as nx

D = "[diretório_transcriptoma]"
DEG = "[diretório_deseq2_results]/DEGs_todos_shrink.tsv"
tid = lambda s: s.str.replace(r"\.p\d+$", "", regex=True)
dg = pd.read_csv(DEG, sep="\t").set_index("gene")

# 1. diagnóstico: proporção up/down na entrada vs na rede
h = pd.read_csv(f"{D}/diamond_agave_canonical_degs_best_hits.tsv", sep="\t", header=None, usecols=[0], names=["t"])
l = dg["log2FoldChange"].reindex(tid(h.t)).values
print(f"Entrada canônica: {len(h)} transcritos | up = {np.nansum(l > 0)} ({100*np.nansum(l > 0)/np.sum(~np.isnan(l)):.1f}%) | down = {np.nansum(l < 0)}")

net = pd.read_csv(f"{D}/rede_canais_permitidos_thr0.7.tsv", sep="\t")
nodes = pd.Series(pd.concat([net.agave_A, net.agave_B]).unique())
l2 = dg["log2FoldChange"].reindex(tid(nodes)).values
print(f"Nós com aresta em 0,7 (antes de qualquer filtro): {len(nodes)} | up = {np.sum(l2 > 0)} ({100*np.sum(l2 > 0)/len(l2):.1f}%)")

# 2. representante por referência
ref = {}
for a, b, ra, rb in zip(net.agave_A, net.agave_B, net.ref_A, net.ref_B):
    ref[a] = ra; ref[b] = rb
t = pd.DataFrame({"node": list(ref), "ref": list(ref.values())})
t["transcript_id"] = tid(t.node)
t["log2FoldChange"] = dg["log2FoldChange"].reindex(t.transcript_id).values
t["padj"] = dg["padj"].reindex(t.transcript_id).values
t["baseMean"] = dg["baseMean"].reindex(t.transcript_id).values
t["abs"] = t.log2FoldChange.abs()
g = t.groupby("ref")
resumo = pd.DataFrame({"n_transcritos": g.size(), "n_up": g.apply(lambda x: (x.log2FoldChange > 0).sum()),
                       "n_down": g.apply(lambda x: (x.log2FoldChange < 0).sum())})
resumo["conflito"] = (resumo.n_up > 0) & (resumo.n_down > 0)
print(f"Referências com isoformas em direções opostas: {int(resumo.conflito.sum())} de {len(resumo)}")
rep = t.sort_values(["padj", "abs"], ascending=[True, False]).drop_duplicates("ref").set_index("ref")
rep = rep.join(resumo)

# 3. rede colapsada (nós = proteínas de referência)
e = net[["ref_A", "ref_B", "score_permitido"]].copy()
e[["ref_A", "ref_B"]] = np.sort(e[["ref_A", "ref_B"]].values, axis=1)
e = e[e.ref_A != e.ref_B].groupby(["ref_A", "ref_B"], as_index=False)["score_permitido"].max()
G = nx.from_pandas_edgelist(e, "ref_A", "ref_B", edge_attr="score_permitido")
for c in [c for c in nx.connected_components(G) if len(c) < 3]:
    G.remove_nodes_from(c)
print(f"Rede colapsada: {G.number_of_nodes()} nós, {G.number_of_edges()} arestas")

m = rep.loc[list(G.nodes())].copy()
m["degree"] = pd.Series(dict(G.degree()))
m["betweenness"] = pd.Series(nx.betweenness_centrality(G))
com = sorted(nx.community.louvain_communities(G, seed=42), key=len, reverse=True)
m["module"] = pd.Series({n: i + 1 for i, c in enumerate(com) for n in c})
m["direction"] = np.where(m.log2FoldChange > 0, "up", "down")
print(f"Up: {(m.direction == 'up').sum()} | Down: {(m.direction == 'down').sum()}")
print("\nMódulos (>= 10 nós):")
print(m.groupby("module").agg(nos=("degree", "size"), up=("direction", lambda s: (s == "up").sum()),
      down=("direction", lambda s: (s == "down").sum()), grau_max=("degree", "max")).query("nos >= 10").to_string())

m.sort_values("degree", ascending=False).to_csv(f"{D}/rede_colapsada_nos.tsv", sep="\t")
nx.to_pandas_edgelist(G).rename(columns={"source": "ref_A", "target": "ref_B"}).to_csv(f"{D}/rede_colapsada_arestas.tsv", sep="\t", index=False)
print("Pronto.")
