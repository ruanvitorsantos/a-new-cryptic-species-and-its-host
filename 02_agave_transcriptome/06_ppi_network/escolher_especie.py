# Protein-protein interaction network analysis (Methods 5.2.6): escolher especie.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import numpy as np, pandas as pd, networkx as nx
D = "[diretório_transcriptoma]"
DEG = "[diretório_deseq2_results]/DEGs_todos_shrink.tsv"
NOME = {"4615": "abacaxi", "906689": "dendrobium", "3702": "arabidopsis", "4577": "milho", "39947": "arroz", "4558": "sorgo", "214687": "banana"}
dg = pd.read_csv(DEG, sep="\t").set_index("gene")["log2FoldChange"]
N = len(open(f"{D}/ref5/ids_canonicos.txt").read().split())
h = pd.read_csv(f"{D}/ref7/hits7_all.tsv", sep="\t", header=None, names=["q", "s", "e", "b", "p", "l"])
h["tax"] = h.s.str.split(".").str[0]
best = h.sort_values("b", ascending=False).drop_duplicates(["q", "tax"])
E = pd.read_csv(f"{D}/consenso7/arestas_consenso_todas.tsv", sep="\t")
E["sp"] = E.referencias.str.split(";")
cons = E[E.K >= 3]
rows = []
for t, n in NOME.items():
    b = best[best.tax == t]
    sub = E[E.sp.apply(lambda l: n in l)]
    G = nx.from_pandas_edgelist(sub, "A", "B")
    for c in [c for c in nx.connected_components(G) if len(c) < 3]:
        G.remove_nodes_from(c)
    l = dg.reindex(list(G.nodes()))
    rec = cons.sp.apply(lambda x: n in x).sum() / len(cons)
    rows.append(dict(especie=n, cobertura_pct=round(100 * b.q.nunique() / N, 1), pident_mediana=round(b.p.median(), 1),
                     nos=G.number_of_nodes(), arestas=G.number_of_edges(), up=int((l > 0).sum()),
                     pct_up=round(100 * (l > 0).mean(), 1), recupera_consenso_pct=round(100 * rec, 1),
                     arestas_com_K3_pct=round(100 * (sub.K >= 3).mean(), 1)))
R = pd.DataFrame(rows).sort_values("recupera_consenso_pct", ascending=False)
print(R.to_string(index=False))
R.to_csv(f"{D}/escolha_especie.tsv", sep="\t", index=False)
