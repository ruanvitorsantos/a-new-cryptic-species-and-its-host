#!/bin/bash
# Protein-level overlap between the two transcriptomes (DIAMOND blastp, >= 90% identity, e-value 1e-5).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


# =====================================================================
# AJUSTE AQUI: Coloque os nomes exatos dos seus arquivos de PROTEÍNAS
MINHA_PROTEINA="Agave_sisalana_transcriptome.fasta.pep"
PROTEINA_UNICAMP="Transcriptome_Codings_Agave_sisalana.fasta.pep"
# =====================================================================

echo "--------------------------------------------------"
echo "1/2: Criando o Banco de Dados de Proteínas (Unicamp)..."
echo "--------------------------------------------------"
diamond makedb --in "$PROTEINA_UNICAMP" -d db_unicamp_prot

echo "--------------------------------------------------"
echo "2/2: Iniciando o alinhamento ULTRA RÁPIDO (blastp)..."
echo "--------------------------------------------------"
diamond blastp \
        -q "$MINHA_PROTEINA" \
        -d db_unicamp_prot \
        -o overlap_blast/resultado_overlap_diamond.txt \
        -p 8 \
        -e 1e-5 \
        --id 90 \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "--------------------------------------------------"
echo "Overlap concluído com sucesso!"
echo "Resultado salvo em: overlap_blast/resultado_overlap_diamond.txt"
echo "--------------------------------------------------"
