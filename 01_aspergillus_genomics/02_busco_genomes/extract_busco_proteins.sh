#!/usr/bin/env bash
# Alternative/complementary step: collect the single-copy BUSCO proteins of every
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

# finished BUSCO run into one directory (one .faa per genome, headers prefixed with the strain ID).
# Run inside the directory that holds the BUSCO result folders.

DESTINO="[diretório_proteinas_busco]/"
mkdir -p "${DESTINO}"

for pasta in */run_aspergillus_odb12/; do
    SAMPLE=$(echo "${pasta}" | cut -d'/' -f1)
    echo "Processando: ${SAMPLE}"
    ARQUIVO_SAIDA="${DESTINO}${SAMPLE}.faa"

    cat "${pasta}busco_sequences/single_copy_busco_sequences/"*faa > "${ARQUIVO_SAIDA}" 2>/dev/null

    if [ -s "${ARQUIVO_SAIDA}" ]; then
        sed -i 's,>,>'"${SAMPLE}"'_,' "${ARQUIVO_SAIDA}"
    else
        echo "Aviso: Nenhuma sequência encontrada para ${SAMPLE}"
        rm -f "${ARQUIVO_SAIDA}"
    fi
done

echo "Concluído! Todas as proteínas foram extraídas para ${DESTINO}"
