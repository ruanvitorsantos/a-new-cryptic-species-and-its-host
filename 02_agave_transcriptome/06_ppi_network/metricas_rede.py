# Protein-protein interaction network analysis (Methods 5.2.6): metricas rede.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import re, numpy as np, pandas as pd, networkx as nx

D = "[diretório_transcriptoma]"
DEG = "[diretório_deseq2_results]/DEGs_todos_shrink.tsv"
net = pd.read_csv(f"{D}/rede_canais_permitidos_thr0.7.tsv", sep="\t")
G = nx.from_pandas_edgelist(net, "agave_A", "agave_B", edge_attr="score_permitido")

# mantém só componentes com >= 3 nós (regra da legenda do Frontiers)
for c in [c for c in nx.connected_components(G) if len(c) < 3]:
    G.remove_nodes_from(c)
print(f"Após remover componentes < 3 nós: {G.number_of_nodes()} nós, {G.number_of_edges()} arestas")

# nó -> proteína de referência
ref = {}
for a, b, ra, rb in zip(net.agave_A, net.agave_B, net.ref_A, net.ref_B):
    ref[a] = ra; ref[b] = rb
nref = pd.Series({n: ref[n] for n in G.nodes()})
por_ref = nref.value_counts()
print(f"Proteínas de referência distintas: {por_ref.size} | referências com >1 transcrito: {(por_ref > 1).sum()} "
      f"| nós em referências compartilhadas: {int(por_ref[por_ref > 1].sum())} | máximo de transcritos por referência: {por_ref.max()}")

# vizinhança idêntica (gêmeos topológicos)
viz = {n: frozenset(G[n]) - {n} for n in G}
assinatura = pd.Series({n: (nref[n], viz[n]) for n in G}).astype(str)
print(f"Nós gêmeos (mesma referência e mesma vizinhança): {int(assinatura.duplicated(keep=False).sum())}")

# métricas e módulos
deg = pd.Series(dict(G.degree()))
bet = pd.Series(nx.betweenness_centrality(G))
com = sorted(nx.community.louvain_communities(G, seed=42), key=len, reverse=True)
mod = {n: i + 1 for i, c in enumerate(com) for n in c}

dg = pd.read_csv(DEG, sep="\t").set_index("gene")
m = pd.DataFrame({"degree": deg, "betweenness": bet, "module": pd.Series(mod), "ref": nref})
m["transcript_id"] = [re.sub(r"\.p\d+$", "", n) for n in m.index]
m["log2FoldChange"] = dg["log2FoldChange"].reindex(m.transcript_id).values
m["padj"] = dg["padj"].reindex(m.transcript_id).values
m["direction"] = np.where(m.log2FoldChange > 0, "up", "down")
print(f"Nós sem DEG na tabela: {int(m.log2FoldChange.isna().sum())}")
print(f"Up: {(m.direction == 'up').sum()} | Down: {(m.direction == 'down').sum()}")
print("\nMódulos (>= 10 nós):")
print(m.groupby("module").agg(nos=("degree", "size"), up=("direction", lambda s: (s == "up").sum()),
                              down=("direction", lambda s: (s == "down").sum()),
                              grau_max=("degree", "max")).query("nos >= 10").to_string())

m.sort_values("degree", ascending=False).to_csv(f"{D}/rede_canais_nos_metricas.tsv", sep="\t")
nx.to_pandas_edgelist(G).rename(columns={"source": "agave_A", "target": "agave_B"}).to_csv(
    f"{D}/rede_canais_arestas_cytoscape.tsv", sep="\t", index=False)
print("Pronto.")
