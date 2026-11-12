#!/usr/bin/env bash
# ANI by alignment (ANIm, MUMmer) with pyANI v0.2.13.1 (Section 4.2.4; Supplementary Table 4).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


# environment
conda create -n pyani_env python=3.9 -y
conda activate pyani_env
pip install pyani==0.2.13.1
# if the first run complains about missing modules:
pip install "numpy<2" "matplotlib<3.6" biopython intervaltree namedlist networkx openpyxl pandas Pillow pyarrow scipy seaborn sqlalchemy tqdm

# ANIm (MUMmer-based; recommended for closely related genomes)
average_nucleotide_identity.py -i genomas -o ANI_analysis/resultados/ANIm -m ANIm -f --workers 12

# For reference, the BLAST-based variant (ANIb) is run with  -m ANIb  (not used in the manuscript).
