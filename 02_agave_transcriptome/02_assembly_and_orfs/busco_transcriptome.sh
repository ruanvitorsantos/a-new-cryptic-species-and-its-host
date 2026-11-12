#!/bin/bash
# BUSCO (transcriptome mode), lineages liliopsida_odb12 and viridiplantae_odb12, for the transcriptome of this study and the external transcriptome.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


# 1. Liliopsida - Transcriptoma (este estudo)
echo "Iniciando: Liliopsida no transcriptoma deste estudo..."
busco -m transcriptome \
      -i Agave_sisalana_transcriptome.fasta \
      -o busco/resultado_busco_liliopsida_transcriptoma \
      -l [diretório_transcriptoma_base]/busco/lineages/liliopsida_odb12 \
      -c 12 \
      -f

# 2. Liliopsida - Transcriptoma Externo (Unicamp)
echo "Iniciando: Liliopsida no transcriptoma da Unicamp..."
busco -m transcriptome \
      -i Transcriptome_Codings_Agave_sisalana.fasta \
      -o busco/resultado_busco_liliopsida_unicamp \
      -l [diretório_transcriptoma_base]/busco/lineages/liliopsida_odb12 \
      -c 12 \
      -f

# 3. Viridiplantae - Transcriptoma (este estudo)
echo "Iniciando: Viridiplantae no transcriptoma deste estudo..."
busco -m transcriptome \
      -i Agave_sisalana_transcriptome.fasta \
      -o busco/resultado_busco_viridiplantae_transcriptoma \
      -l [diretório_transcriptoma_base]/busco/lineages/viridiplantae_odb12 \
      -c 12 \
      -f

# 4. Viridiplantae - Transcriptoma Externo (Unicamp)
echo "Iniciando: Viridiplantae no transcriptoma da Unicamp..."
busco -m transcriptome \
      -i Transcriptome_Codings_Agave_sisalana.fasta \
      -o busco/resultado_busco_viridiplantae_unicamp \
      -l [diretório_transcriptoma_base]/busco/lineages/viridiplantae_odb12 \
      -c 12 \
      -f

echo "Todos os testes do BUSCO foram concluídos com sucesso!"
