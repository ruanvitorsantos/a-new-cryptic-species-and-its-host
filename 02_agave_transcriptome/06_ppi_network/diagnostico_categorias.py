# Protein-protein interaction network analysis (Methods 5.2.6): diagnostico categorias.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import re, numpy as np, pandas as pd
from functools import lru_cache

D = "[diretório_transcriptoma]"
IPS = re.search(r'^IPS = "(.*)"', open(f"{D}/subrede_visivel.py").read(), re.M).group(1)
DEG = "[diretório_deseq2_results]/DEGs_todos_shrink.tsv"
PEP = f"{D}/protein/Agave_sisalana_transcriptoma_sem_rRNA.fasta.transdecoder.pep"

par, nome, cur = {}, {}, None
for l in open(f"{D}/go-basic.obo"):
    l = l.strip()
    if l.startswith("["): cur = None
    if l == "[Term]": cur = "?"
    elif l.startswith("id: GO:") and cur == "?": cur = l[4:]; par[cur] = set()
    elif cur and cur != "?" and l.startswith("name: "): nome[cur] = l[6:]
    elif cur and cur != "?" and l.startswith("is_a: GO:"): par[cur].add(l[6:16])
    elif cur and cur != "?" and l.startswith("relationship: part_of GO:"): par[cur].add(l[22:32])
@lru_cache(None)
def anc(g): return frozenset({g}).union(*[anc(p) for p in par.get(g, ())]) if g in par else frozenset({g})
print("Teste do parser:")
for g in ("GO:0006950", "GO:0009408", "GO:0006457", "GO:0006412"):
    print(f"  {g} ({nome.get(g)}): {len(anc(g))} ancestrais | contém GO:0050896 = {'GO:0050896' in anc(g)}")

CAT = {"estresse/estímulo": ["GO:0050896"], "proteínas (dobramento/degradação)": ["GO:0006457", "GO:0006508", "GO:0043161", "GO:0000502"],
       "tradução/ribossomo": ["GO:0006412", "GO:0042254", "GO:0005840"], "carboidratos/crescimento": ["GO:0005975", "GO:0040007"],
       "DNA": ["GO:0006260", "GO:0006281"]}
ips = pd.read_csv(IPS, sep="\t", header=None, usecols=[0, 13], names=["p", "go"], dtype=str).dropna()
ips["go"] = ips.go.str.findall(r"GO:\d{7}")
ips["t"] = ips.p.str.replace(r"\.p\d+$", "", regex=True)
tg = ips.explode("go").dropna(subset=["go"]).groupby("t")["go"].agg(set)

orf = set()
for l in open(PEP):
    if l[0] == ">": orf.add(re.sub(r"\.p\d+$", "", l[1:].split()[0]))
dg = pd.read_csv(DEG, sep="\t").set_index("gene")
d = dg[(dg.padj < 0.05) & (dg.log2FoldChange.abs() > 1) & dg.index.isin(orf)].copy()
d["dir"] = np.where(d.log2FoldChange > 0, "up", "down")
d["tem_go"] = d.index.isin(tg.index)
def cats(t):
    A = set()
    for g in tg.get(t, set()): A |= anc(g)
    return {c for c, rs in CAT.items() if A & set(rs)}
d["cats"] = [cats(t) if h else set() for t, h in zip(d.index, d.tem_go)]
out = {}
for dr, x in d.groupby("dir"):
    out[dr] = {"transcritos com ORF": len(x), "% com algum GO": 100 * x.tem_go.mean(),
               **{f"% {c}": 100 * np.mean([c in s for s in x.cats]) for c in CAT}}
print("\nTodos os DEGs com ORF:"); print(pd.DataFrame(out).round(1).to_string())

for dr in ("up", "down"):
    ts = d.index[(d.dir == dr) & d.tem_go]
    cnt = pd.Series([g for t in ts for g in tg[t]]).value_counts().head(12)
    print(f"\nTermos GO mais frequentes ({dr}):")
    for g, n in cnt.items(): print(f"  {n:5d}  {g}  {nome.get(g, '?')}")

s = pd.read_csv(f"{D}/subrede_150_nos.tsv", sep="\t")
s["tem_go"] = s.transcript_id.isin(tg.index)
print("\nNós da sub-rede com algum GO (por direção):"); print(s.groupby("direction").tem_go.agg(["size", "mean"]).round(2).to_string())
