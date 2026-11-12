# Protein-protein interaction network analysis (Methods 5.2.6): recuperar up.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import os, re, numpy as np, pandas as pd
D = "[diretório_transcriptoma]"
Q = "[diretório_salmon_quants]"
DEG = "[diretório_deseq2_results]/DEGs_todos_shrink.tsv"
MINLEN = int(os.environ.get("MINLEN", "30"))
SD = ["ASD2_concatenado_S28_S17_S6", "ASD3_concatenado_S29_S18_S7", "ASD5_concatenado_S30_S19_S8",
      "ASD6_concatenado_S31_S20_S9", "ASD7_concatenado_S32_S21_S10", "ASD8_concatenado_S33_S22_S11"]
SS = ["ASS2_concatenado_S23_S12_S1", "ASS5_concatenado_S24_S13_S2", "ASS7_concatenado_S25_S14_S3",
      "ASS8_concatenado_S26_S15_S4", "ASS9_concatenado_S27_S16_S5"]
gene = lambda s: re.sub(r"_i\d+(\.p\d+)?$", "", s)
dg = pd.read_csv(DEG, sep="\t").set_index("gene")
deg = dg[(dg.padj < 0.05) & (dg.log2FoldChange.abs() > 1)].copy()
deg["dir"] = np.where(deg.log2FoldChange > 0, "up", "down")
cont = lambda ids: deg["dir"].reindex(list(ids)).value_counts().reindex(["up", "down"]).fillna(0).astype(int)
h = pd.read_csv(f"{D}/ref7/semorf_todos_blastx.tsv", sep="\t", header=None, names=["q", "s", "e", "b", "p", "l"])
print("Sem ORF com hit, sem filtro de comprimento:", cont(h.q.unique()).to_dict())
h = h[h.l >= MINLEN]
cand = h.q.unique()
tp = pd.DataFrame({a: pd.read_csv(f"{Q}/clean_{a}_matched/quant.sf", sep="\t", usecols=["Name", "TPM"]).set_index("Name")["TPM"].reindex(cand).fillna(0) for a in SD + SS})
ok = (tp[SD].mean(axis=1) >= 1) | (tp[SS].mean(axis=1) >= 1)
ct = [re.sub(r"\.p\d+$", "", x) for x in open(f"{D}/ref5/ids_canonicos.txt").read().split()]
cg = {gene(x) for x in ct}
best = h.sort_values("b", ascending=False).drop_duplicates("q").set_index("q")
best["tpm"] = tp.sum(axis=1).reindex(best.index); best["gene"] = [gene(x) for x in best.index]
best["tpm_ok"] = best.index.isin(tp.index[ok]); best["gene_novo"] = ~best.gene.isin(cg)
sel = best[best.tpm_ok & best.gene_novo].sort_values(["b", "tpm"], ascending=False).drop_duplicates("gene")
t = pd.DataFrame({f"com hit (>= {MINLEN} aa)": cont(best.index), "+ TPM >= 1": cont(best.index[best.tpm_ok]),
                  "+ gene fora do conjunto canônico": cont(best.index[best.tpm_ok & best.gene_novo]), "recuperados (1 por gene)": cont(sel.index)})
print(t.T.to_string())
c0 = cont(ct); c1 = c0 + cont(sel.index)
print("Conjunto canônico antes:", c0.to_dict(), f"({100 * c0['up'] / c0.sum():.1f}% up) | depois:", c1.to_dict(), f"({100 * c1['up'] / c1.sum():.1f}% up)")
print("Espécie do melhor hit dos recuperados:", sel.s.str.split(".").str[0].value_counts().to_dict())
sel.assign(dir=deg["dir"].reindex(sel.index).values, log2FC=dg.log2FoldChange.reindex(sel.index).values).to_csv(f"{D}/ref7/recuperados_info.tsv", sep="\t")
orig = pd.read_csv(f"{D}/ref7/hits7_all.tsv", sep="\t", header=None, names=["q", "s", "e", "b", "p", "l"])
pd.concat([orig, h.loc[h.q.isin(sel.index), ["q", "s", "e", "b", "p", "l"]]]).to_csv(f"{D}/ref7/hits7_all_plus.tsv", sep="\t", header=False, index=False)
print("Gravados: ref7/recuperados_info.tsv e ref7/hits7_all_plus.tsv")
