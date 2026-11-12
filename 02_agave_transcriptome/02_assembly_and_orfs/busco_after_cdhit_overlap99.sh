#!/bin/bash
# BUSCO on the 99% CD-HIT overlap set.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


# 1. Liliopsida - Transcriptoma Externo (Unicamp)
echo "Iniciando: Liliopsida.."
busco -m transcriptome \
      -i transcriptoma_overlap_99 \
      -o busco_overlap_99/resultado_busco_liliopsida \
      -l [diretório_transcriptoma_base]/busco/lineages/liliopsida_odb12 \
      -c 12 \
      -f

# 2. Viridiplantae - Transcriptoma (este estudo)
echo "Iniciando: Viridiplantae..."
busco -m transcriptome \
      -i transcriptoma_overlap_99 \
      -o busco_overlap_99/resultado_busco_viridiplantae \
      -l [diretório_transcriptoma_base]/busco/lineages/viridiplantae_odb12 \
      -c 12 \
      -f

echo "Todos os testes do BUSCO foram concluídos com sucesso!"
