# Protein-protein interaction network analysis (Methods 5.2.6): genes e funcoes.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import gzip, numpy as np, pandas as pd
from scipy.stats import hypergeom
D = "[diretório_transcriptoma]"
TAX = {"arabidopsis": "3702", "abacaxi": "4615", "arroz": "39947", "milho": "4577", "sorgo": "4558", "banana": "214687", "dendrobium": "906689"}
pasta = lambda t: D if t in ("4615", "906689") else f"{D}/string_db5"
a = pd.read_csv(f"{D}/rede_visual_nos.tsv", sep="\t").set_index("transcript_id")
b = pd.read_csv(f"{D}/rede_consenso_K3_nos.tsv", sep="\t", index_col=0)
full = pd.read_csv(f"{D}/consenso7/nos_consenso_K3.tsv", sep="\t", index_col=0)
edges = pd.read_csv(f"{D}/consenso7/arestas_consenso_K3_todas.tsv", sep="\t") if False else pd.read_csv(f"{D}/consenso7/arestas_consenso_todas.tsv", sep="\t")
edges = edges[(edges.K >= 3) & edges.A.isin(full.index) & edges.B.isin(full.index)]

need = {r for s in TAX for r in full[s].dropna()}
names, terms = {}, {}
for t in TAX.values():
    for l in gzip.open(f"{D}/string_db/{t}.protein.info.v12.0.txt.gz" if t in ("4615", "906689") else f"{D}/string_db5/{t}.protein.info.v12.0.txt.gz", "rt"):
        p = l.rstrip("\n").split("\t")
        if p[0] in need: names[p[0]] = p[3] if len(p) > 3 else ""
    with gzip.open(f"{pasta(t)}/{t}.protein.enrichment.terms.v12.0.txt.gz", "rt") as f:
        for l in f:
            if "Biological Process" not in l: continue
            p = l.rstrip("\n").split("\t")
            if p[0] in need: terms[p[0]] = terms.get(p[0], set()) | {(p[2], p[3])}
def anotacao(n):
    for s in TAX:
        r = full.at[n, s]
        if pd.notna(r) and names.get(r, "") and not names[r].lower().startswith("uncharacterized"): return f"{s}: {names[r]}"
    return "sem anotação informativa"

# 1. genes up compartilhados e exclusivos
ua, ub = set(a.index[a.direction == "up"]), set(b.index[b.direction == "up"])
rows = []
for grp, S, fonte in (("compartilhado", ua & ub, b), ("só A (abacaxi+dendrobium)", ua - ub, a), ("só B (consenso)", ub - ua, b)):
    for n in sorted(S, key=lambda x: -fonte.at[x, "log2FoldChange"]):
        rows.append(dict(grupo=grp, transcrito=n, log2FC=round(fonte.at[n, "log2FoldChange"], 2), rotulo=fonte.at[n, "rotulo"],
                         anotacao=anotacao(n) if n in full.index else "fora da rede de consenso", categorias=fonte.at[n, "cats"]))
T = pd.DataFrame(rows); T.to_csv(f"{D}/genes_up_A_vs_B.tsv", sep="\t", index=False)
print(T.grupo.value_counts().to_string()); print(T.to_string(index=False, max_colwidth=60))

# 2. estresse só com GO das referências que sustentam as arestas
ST = {"GO:0006950", "GO:0006952"}
par, cur = {}, None
for l in open(f"{D}/go-basic.obo"):
    l = l.strip()
    if l.startswith("["): cur = None
    if l == "[Term]": cur = "?"
    elif l.startswith("id: GO:") and cur == "?": cur = l[4:]; par[cur] = set()
    elif cur and cur != "?" and l.startswith("is_a: GO:"): par[cur].add(l[6:16])
    elif cur and cur != "?" and l.startswith("relationship: part_of GO:"): par[cur].add(l[22:32])
from functools import lru_cache
@lru_cache(None)
def anc(g): return frozenset({g}).union(*[anc(p) for p in par.get(g, ())]) if g in par else frozenset({g})
def tem_estresse(n, so_sup):
    sup = set()
    if so_sup:
        for refs in edges[(edges.A == n) | (edges.B == n)].referencias: sup |= set(refs.split(";"))
    A = set()
    for s in TAX:
        r = full.at[n, s]
        if pd.notna(r) and (not so_sup or s in sup):
            for g, _ in terms.get(r, ()): A |= anc(g)
    return bool(A & ST)
sh = sorted(set(a.index) & set(b.index))
cA = a.loc[sh, "cats"].fillna("").str.contains("Resposta a estresse", regex=False)
cU = pd.Series([tem_estresse(n, False) for n in sh], index=sh); cS = pd.Series([tem_estresse(n, True) for n in sh], index=sh)
print(f"\nNós em comum (n={len(sh)}): estresse A={int(cA.sum())} | B união das 7 refs={int(cU.sum())} | B só refs que sustentam as arestas={int(cS.sum())}")

# 3. enriquecimento por módulo (B, K>=3), universo = nós da rede K>=3
GOT = {}
for n in full.index:
    s = set()
    for sp in TAX:
        r = full.at[n, sp]
        if pd.notna(r): s |= {g for g, _ in terms.get(r, ())}
    GOT[n] = s
nomego = {g: d for v in terms.values() for g, d in v}
U = list(full.index); NU = len(U)
cont = {}
for n in U:
    for g in GOT[n]: cont[g] = cont.get(g, 0) + 1
out = []
for mo, grp in full.groupby("module"):
    if len(grp) < 15: continue
    cm = {}
    for n in grp.index:
        for g in GOT[n]: cm[g] = cm.get(g, 0) + 1
    for g, k in cm.items():
        if k < 3 or cont[g] > 0.2 * NU: continue
        out.append(dict(modulo=mo, nos_modulo=len(grp), go=g, termo=nomego.get(g, ""), k=k, K_universo=cont[g], p=hypergeom.sf(k - 1, NU, cont[g], len(grp))))
E = pd.DataFrame(out)
E["fdr"] = np.minimum.accumulate((E.p.sort_values(ascending=False) * len(E) / (len(E) - np.arange(len(E)))).values)[::-1] if False else None
o = E.sort_values("p").reset_index(drop=True); o["fdr"] = np.minimum.accumulate((o.p * len(o) / (np.arange(len(o)) + 1))[::-1])[::-1].clip(upper=1)
o = o[o.fdr < 0.05]; o.to_csv(f"{D}/enriquecimento_modulos_K3.tsv", sep="\t", index=False)
ud = full.groupby("module").direction.value_counts().unstack(fill_value=0)
print("\nMódulos (nós, up, down):"); print(ud.assign(n=ud.sum(axis=1)).query("n >= 15").to_string())
print("\nTermos enriquecidos (FDR < 0,05) por módulo, top 5:")
for mo, g in o.groupby("modulo"): print(f"  módulo {mo}:", "; ".join(f"{r.termo} ({r.k}/{r.nos_modulo})" for r in g.head(5).itertuples()))
