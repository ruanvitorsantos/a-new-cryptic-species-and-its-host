"""Summarise % identity and mismatches of a BLAST outfmt-6 table (Agave transcriptome vs Asparagus officinalis).
Usage: python3 summarize_blast_identity.py resultado_agave_vs_aspargo.tsv
"""
import sys

arquivo = sys.argv[1] if len(sys.argv) > 1 else "resultado_agave_vs_aspargo.tsv"
identidades, mismatches = [], []

try:
    with open(arquivo) as f:
        for linha in f:
            colunas = linha.strip().split("\t")
            if len(colunas) >= 5:
                identidades.append(float(colunas[2]))   # column 3: % identity
                mismatches.append(int(colunas[4]))      # column 5: absolute mismatches

    if identidades:
        print("\n" + "=" * 55)
        print("TOLERANCE ANALYSIS (AGAVE VS ASPARAGUS)")
        print("=" * 55)
        print(f"Highest identity (most conserved):  {max(identidades):.2f}%")
        print(f"Lowest identity (least conserved):  {min(identidades):.2f}%")
        print(f"Mean identity:                      {sum(identidades)/len(identidades):.2f}%")
        print("-" * 55)
        print(f"Maximum mismatches found:           {max(mismatches)}")
        print(f"Mean mismatches per alignment:      {sum(mismatches)/len(mismatches):.1f}")
        print("=" * 55 + "\n")
    else:
        print("The file is empty: no BLAST hit was found.")
except FileNotFoundError:
    print(f"Error: {arquivo} not found. Make sure BLAST has finished.")
