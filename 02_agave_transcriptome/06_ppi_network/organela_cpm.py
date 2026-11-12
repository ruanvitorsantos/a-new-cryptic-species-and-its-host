# Protein-protein interaction network analysis (Methods 5.2.6): organela cpm.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import re, pandas as pd, numpy as np
from scipy.stats import mannwhitneyu
D = "[diretório_transcriptoma]"
Q = "[diretório_salmon_quants]"
AM = {"ASD2": "ASD2_concatenado_S28_S17_S6", "ASD3": "ASD3_concatenado_S29_S18_S7", "ASD5": "ASD5_concatenado_S30_S19_S8",
      "ASD6": "ASD6_concatenado_S31_S20_S9", "ASD7": "ASD7_concatenado_S32_S21_S10", "ASD8": "ASD8_concatenado_S33_S22_S11",
      "ASS2": "ASS2_concatenado_S23_S12_S1", "ASS5": "ASS5_concatenado_S24_S13_S2", "ASS7": "ASS7_concatenado_S25_S14_S3",
      "ASS8": "ASS8_concatenado_S26_S15_S4", "ASS9": "ASS9_concatenado_S27_S16_S5"}
PAD = r"^(atp|nad|ndh|rps|rpl|psa|psb|pet|rbc|cox|cob|ccm|mat|rpo|ycf)"
full = pd.read_csv(f"{D}/consenso7/nos_consenso_K3.tsv", sep="\t", index_col=0)
inf = pd.read_csv(f"{D}/string_db5/3702.protein.info.v12.0.txt.gz", sep="\t", header=0, usecols=[0, 1], names=["id", "sym"]).set_index("id")["sym"]
ara = full["arabidopsis"].dropna()
org = pd.Series([bool(re.match(PAD, str(s))) for s in inf.reindex(ara).values], index=ara.index)
ids = list(org[org].index)
print("Transcritos de organela (critério de símbolo):", len(ids))
cpm, frac = {}, {}
for a, f in AM.items():
    q = pd.read_csv(f"{Q}/clean_{f}_matched/quant.sf", sep="\t", usecols=["Name", "NumReads"]).set_index("Name")["NumReads"]
    tot = q.sum()
    x = q.reindex(ids).fillna(0)
    cpm[a] = x / tot * 1e6
    frac[a] = 100 * x.sum() / tot
F = pd.Series(frac)
print("\n% das leituras em transcritos de organela:")
print(F.round(3).to_string())
d = F[[k for k in F.index if k.startswith("ASD")]]
s = F[[k for k in F.index if k.startswith("ASS")]]
print(f"\nMediana doente = {d.median():.3f}% | sadia = {s.median():.3f}% | Mann-Whitney p = {mannwhitneyu(d, s)[1]:.3g}")
C = pd.DataFrame(cpm)
cd = C[[c for c in C if c.startswith("ASD")]].median(axis=1)
cs = C[[c for c in C if c.startswith("ASS")]].median(axis=1)
r = np.log2((cd + 1) / (cs + 1))
print("\nlog2(CPM doente/sadia), transcritos de organela:")
print(f"  mediana = {r.median():.2f} | n com log2 > 1: {int((r > 1).sum())} | n com log2 < -1: {int((r < -1).sum())}")
