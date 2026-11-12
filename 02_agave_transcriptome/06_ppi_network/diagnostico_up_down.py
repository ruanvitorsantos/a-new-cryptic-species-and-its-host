# Protein-protein interaction network analysis (Methods 5.2.6): diagnostico up down.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import numpy as np, pandas as pd

D = "[diretório_transcriptoma]"
DEG = "[diretório_deseq2_results]/DEGs_todos_shrink.tsv"
OK = ["neighborhood", "fusion", "cooccurence", "coexpression", "experimental"]
ALL = OK + ["database", "textmining"]; P = 0.041

def comb(df, cols):
    prod = np.ones(len(df))
    for c in cols:
        prod *= 1 - np.clip((df[c].values - P) / (1 - P), 0, None)
    s = 1 - prod
    return s + P * (1 - s)

h = pd.read_csv(f"{D}/diamond_agave_canonical_degs_best_hits.tsv", sep="\t", header=None, usecols=[0, 1], names=["t", "ref"])
dg = pd.read_csv(DEG, sep="\t").set_index("gene")["log2FoldChange"]
h["lfc"] = dg.reindex(h.t.str.replace(r"\.p\d+$", "", regex=True)).values
h = h.dropna(subset=["lfc"]); h["dir"] = np.where(h.lfc > 0, "up", "down")
refs = set(h.ref); print(f"Transcritos: {len(h)} | referências: {len(refs)}")

res = {}
for t in (4615, 906689):
    print(f"Lendo {t}...")
    for ch in pd.read_csv(f"{D}/string_db/{t}.protein.links.detailed.v12.0.txt.gz", sep=" ", chunksize=2_000_000):
        ch = ch[ch.protein1.isin(refs)].copy()
        if ch.empty: continue
        ch[ALL] = ch[ALL] / 1000
        ch["ok"] = comb(ch, OK)
        ch["ok_set"] = np.where(ch.protein2.isin(refs), ch.ok, 0)
        ch["tem"] = 1
        g = ch.groupby("protein1").agg(tem=("tem", "max"), ok=("ok", "max"), ok_set=("ok_set", "max"))
        for p, r in g.iterrows():
            o = res.get(p, (0, 0, 0)); res[p] = (max(o[0], r.tem), max(o[1], r.ok), max(o[2], r.ok_set))
s = pd.DataFrame(res, index=["tem_parceiro", "ok_max", "ok_max_no_conjunto"]).T
h = h.join(s, on="ref").fillna({"tem_parceiro": 0, "ok_max": 0, "ok_max_no_conjunto": 0})
out = h.groupby("dir").apply(lambda x: pd.Series({
    "transcritos": len(x),
    "% com qualquer parceiro no STRING": 100 * x.tem_parceiro.mean(),
    "% com parceiro, canais permitidos >= 0.4": 100 * (x.ok_max >= 0.4).mean(),
    "% com parceiro, canais permitidos >= 0.7": 100 * (x.ok_max >= 0.7).mean(),
    "% com parceiro >= 0.7 DENTRO do conjunto de DEGs": 100 * (x.ok_max_no_conjunto >= 0.7).mean()}))
print(out.T.round(1).to_string())
