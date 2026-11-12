# Protein-protein interaction network analysis (Methods 5.2.6): enriquecimento modulos abacaxi.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import csv, gzip, os
from math import comb
D = "[diretório_projeto_sisal]"
nos = {}
for r in csv.DictReader(open(f"[diretório_transcriptoma]/rede_abacaxi/nos_completa.tsv", encoding="utf-8"), delimiter="\t"):
    nos[r["id"]] = (r["module"], r["direction"])
fn = "4615.protein.enrichment.terms.v12.0.txt.gz"
f = next(os.path.join(r, fn) for r, _, fs in os.walk(D) if fn in fs)
print("termos:", f)
terms, nome = {}, {}
with gzip.open(f, "rt") as fh:
    for l in fh:
        if "Biological Process" not in l: continue
        p = l.rstrip("\n").split("\t")
        if p[0] in nos:
            terms.setdefault(p[0], set()).add(p[2]); nome[p[2]] = p[3]
NU = len(nos); cont = {}
for n in nos:
    for g in terms.get(n, ()): cont[g] = cont.get(g, 0) + 1
def sf(k, K, n):  # P(X >= k), hipergeometrica exata
    tot = comb(NU, n)
    return sum(comb(K, i) * comb(NU - K, n - i) for i in range(k, min(K, n) + 1)) / tot
mods = {}
for n, (m, d) in nos.items(): mods.setdefault(m, []).append(n)
out = []
for mo, grp in mods.items():
    if len(grp) < 10: continue
    cm = {}
    for n in grp:
        for g in terms.get(n, ()): cm[g] = cm.get(g, 0) + 1
    for g, k in cm.items():
        if k < 3 or cont[g] > 0.2 * NU: continue
        out.append([mo, len(grp), g, nome[g], k, cont[g], sf(k, cont[g], len(grp))])
out.sort(key=lambda r: r[6]); N = len(out)
fdr = [0.0] * N; run = 1.0
for i in range(N - 1, -1, -1):
    run = min(run, out[i][6] * N / (i + 1)); fdr[i] = run
S = [r + [fdr[i]] for i, r in enumerate(out) if fdr[i] < 0.05]
os.makedirs("[diretório_tabelas_suplementares]", exist_ok=True)
with open("[diretório_tabelas_suplementares]/S8_module_enrichment_full.tsv", "w") as o:
    o.write("module\tmodule_nodes\tGO_ID\tterm\tk\tK\tp\tFDR\n")
    for r in S: o.write("\t".join(map(str, r)) + "\n")
print("universo", NU, "| testes", N, "| significativos", len(S))
for mo, grp in sorted(mods.items(), key=lambda x: -len(x[1])):
    if len(grp) < 10: continue
    up = sum(nos[n][1] == "up" for n in grp)
    s = sorted([r for r in S if r[0] == mo], key=lambda r: r[7])
    top = " ; ".join(f"{r[3]} ({r[4]}/{r[5]}, {r[7]:.2e})" for r in s[:3])
    print("|".join(map(str, [mo, len(grp), up, len(grp) - up, len(s), top])))
