#!/usr/bin/env bash
# BUSCO on transcriptomes (this study vs. an external A. sisalana transcriptome), lineages liliopsida_odb12 and viridiplantae_odb12.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


# (once) busco --download liliopsida_odb12 ; busco --download viridiplantae_odb12

for LINHAGEM in liliopsida viridiplantae
do
  echo "Iniciando: ${LINHAGEM} no transcriptoma deste estudo..."
  busco -m transcriptome \
        -i Agave_sisalana_transcriptome.fasta \
        -o busco/resultado_busco_${LINHAGEM}_transcriptoma \
        -l [diretório_transcriptoma]/busco/lineages/${LINHAGEM}_odb12 \
        -c 3 -f

  echo "Iniciando: ${LINHAGEM} no transcriptoma externo..."
  busco -m transcriptome \
        -i Transcriptome_Codings_Agave_sisalana.fasta \
        -o busco/resultado_busco_${LINHAGEM}_externo \
        -l [diretório_transcriptoma]/busco/lineages/${LINHAGEM}_odb12 \
        -c 3 -f
done
echo "Todos os testes do BUSCO foram concluídos."
