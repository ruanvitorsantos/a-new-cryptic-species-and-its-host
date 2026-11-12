# Protein-protein interaction network analysis (Methods 5.2.6): orf up down.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

import re, numpy as np, pandas as pd

D = "[diretório_transcriptoma]"
P = f"{D}/protein/Agave_sisalana_transcriptoma_sem_rRNA.fasta"
DEG = "[diretório_deseq2_results]/DEGs_todos_shrink.tsv"
Q = "[diretório_salmon_quants]/clean_ASD2_concatenado_S28_S17_S6_matched/quant.sf"

def ids(path):
    s = set()
    with open(path) as f:
        for l in f:
            if l[0] == ">": s.add(re.sub(r"\.p\d+$", "", l[1:].split()[0]))
    return s

final = ids(f"{P}.transdecoder.pep")
longest = ids(f"{P}.transdecoder_dir/longest_orfs.pep")
dg = pd.read_csv(DEG, sep="\t").set_index("gene")
deg = dg[(dg.padj < 0.05) & (dg.log2FoldChange.abs() > 1)].copy()
deg["dir"] = np.where(deg.log2FoldChange > 0, "up", "down")
deg["comprimento"] = pd.read_csv(Q, sep="\t", usecols=["Name", "Length"]).set_index("Name")["Length"].reindex(deg.index).values
deg["estado"] = np.where(deg.index.isin(final), "1_ORF final",
                np.where(deg.index.isin(longest), "2_so ORF candidato (descartado)", "3_sem ORF candidato"))
print(f"\nTranscritos com ORF final: {len(final)} | com ORF candidato: {len(longest)}")
print((deg.groupby("dir")["estado"].value_counts(normalize=True) * 100).round(1).to_string())
print(deg.groupby(["dir", "estado"]).agg(n=("comprimento", "size"), mediana_comprimento=("comprimento", "median"),
      mediana_baseMean=("baseMean", "median")).round(1).to_string())
