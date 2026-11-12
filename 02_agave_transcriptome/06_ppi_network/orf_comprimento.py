# Protein-protein interaction network analysis (Methods 5.2.6): orf comprimento.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import re, numpy as np, pandas as pd
D = "[diretório_transcriptoma]"
P = f"{D}/protein/Agave_sisalana_transcriptoma_sem_rRNA.fasta"
DEG = "[diretório_deseq2_results]/DEGs_todos_shrink.tsv"
Q = "[diretório_salmon_quants]/clean_ASD2_concatenado_S28_S17_S6_matched/quant.sf"
final = set()
for l in open(f"{P}.transdecoder.pep"):
    if l[0] == ">":
        final.add(re.sub(r"\.p\d+$", "", l[1:].split()[0]))
dg = pd.read_csv(DEG, sep="\t").set_index("gene")
deg = dg[(dg.padj < 0.05) & (dg.log2FoldChange.abs() > 1)].copy()
deg["dir"] = np.where(deg.log2FoldChange > 0, "up", "down")
deg["len"] = pd.read_csv(Q, sep="\t", usecols=["Name", "Length"]).set_index("Name")["Length"].reindex(deg.index).values
deg["orf"] = deg.index.isin(final)
print("Comprimento do transcrito (nt), mediana: ", deg.groupby("dir")["len"].median().to_dict())
print("% com menos de 400 nt:", (deg["len"] < 400).groupby(deg.dir).mean().mul(100).round(1).to_dict())
deg["faixa"] = pd.cut(deg["len"], [0, 300, 400, 500, 700, 1000, 2000, 1e9], labels=["<300", "300-400", "400-500", "500-700", "700-1000", "1000-2000", ">2000"])
taxa = deg.groupby(["faixa", "dir"], observed=True).orf.mean().mul(100).round(1).unstack()
n = deg.groupby(["faixa", "dir"], observed=True).size().unstack()
print("\n% com ORF final, por faixa de comprimento:")
print(taxa.join(n, lsuffix="_%", rsuffix="_n").to_string())
obs = deg[deg.dir == "up"].orf.mean() * 100
esp = (taxa["down"] * n["up"]).sum() / n["up"].sum()
print(f"\nUp com ORF: observado = {obs:.1f}% | esperado se os up tivessem a taxa dos down no mesmo tamanho = {esp:.1f}%")
lens = [int(re.search(r"len:(\d+)", l).group(1)) for l in open(f"{P}.transdecoder_dir/longest_orfs.pep") if l[0] == ">"]
print(f"\nlongest_orfs.pep: {len(lens)} ORFs | menor = {min(lens)} aa | 1o percentil = {np.percentile(lens, 1):.0f} aa | mediana = {np.median(lens):.0f} aa")
