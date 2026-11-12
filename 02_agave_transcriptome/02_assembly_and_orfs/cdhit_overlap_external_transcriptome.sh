#!/bin/bash
# CD-HIT-EST-2D overlap between an external A. sisalana transcriptome and the transcriptome of this study at 90, 95 and 99% identity.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


# Caminho dos arquivos de entrada
DB_UNICAMP="Transcriptome_Codings_Agave_sisalana.fasta"
MINHA_QUERY="Agave_sisalana_transcriptome.fasta"

echo "=========================================================="
echo "1/3: Iniciando Overlap com 90% de Identidade (-c 0.90)"
echo "=========================================================="
cd-hit-est-2d -i "$DB_UNICAMP" \
              -i2 "$MINHA_QUERY" \
              -o cd-hit/transcriptoma_overlap_90 \
              -c 0.9 \
              -n 7 \
              -d 0 \
              -M 0 \
              -T 16

echo "=========================================================="
echo "2/3: Iniciando Overlap com 95% de Identidade (-c 0.95)"
echo "=========================================================="
cd-hit-est-2d -i "$DB_UNICAMP" \
              -i2 "$MINHA_QUERY" \
              -o cd-hit/transcriptoma_overlap_95 \
              -c 0.95 \
              -n 10 \
              -d 0 \
              -M 0 \
              -T 16

echo "=========================================================="
echo "3/3: Iniciando Overlap com 99% de Identidade (-c 0.99)"
echo "=========================================================="
cd-hit-est-2d -i "$DB_UNICAMP" \
              -i2 "$MINHA_QUERY" \
              -o cd-hit/transcriptoma_overlap_99 \
              -c 0.99 \
              -n 10 \
              -d 0 \
              -M 0 \
              -T 16

echo "=========================================================="
echo "Todas as análises do CD-HIT foram concluídas!"
echo "=========================================================="
