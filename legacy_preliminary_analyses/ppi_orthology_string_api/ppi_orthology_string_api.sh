#!/usr/bin/env bash
# =============================================================================
# PIPELINE: Construção de rede de interação proteica (PPI) por ortologia
# para DEGs de transcriptoma de novo (Agave sisalana), via STRING API
#
# Baseado na metodologia de:
# Carvalho LM et al. (2023) "Analysis of protein-protein interaction and
# weighted co-expression networks revealed key modules and genes in
# multiple organs of Agave sisalana." Front. Chem. Eng. 5:1175235.
# DOI: 10.3389/fceng.2023.1175235
#
# Replace the [diretório_*] placeholders with the paths of your environment (see docs/PATHS.md).
# =============================================================================
set -euo pipefail

# --- Project directories --------------------------------------
PROJECT_DIR="[diretório_projeto_sisal]"
TRANSCRIPTOME_DIR="[diretório_transcriptoma]"        # saídas do TransDecoder
DEG_DIR="[diretório_analise_deg]"                    # resultados do DESeq2
ORTHO_DIR="[diretório_rede_ortologia]"                # tudo relacionado a OrthoFinder + STRING

PEP_FILE="${TRANSCRIPTOME_DIR}/transcriptoma_clean.pep"           # proteoma predito (TransDecoder)
DEG_TABLE="${DEG_DIR}/DEGs_resultado_completo.tsv"                 # tabela completa DESeq2 (não filtrada)

mkdir -p "${ORTHO_DIR}"

# =============================================================================
# ETAPA 1 — Verificar se os IDs de proteína são reconhecidos direto pelo STRING
# =============================================================================
# Objetivo: confirmar formato dos headers do proteoma predito.
grep ">" "${PEP_FILE}" | head -n 5
# Resultado esperado: IDs no formato "NOME_TRANSCRITO.p1" (TransDecoder) —
# NÃO reconhecidos diretamente pelo STRING. É necessário mapear por ortologia.

# =============================================================================
# ETAPA 2 — Verificar quais organismos de referência estão disponíveis no STRING
# =============================================================================
# O orientador pediu para checar a espécie mais próxima disponível em:
# https://string-db.org/organism_overview.html

curl -s "https://stringdb-downloads.org/download/species.v12.0.txt" \
    -o /tmp/string_species_v12.txt
wc -l /tmp/string_species_v12.txt

# Checar se a espécie esperada (ex: espécie do mesmo gênero/família) está listada
grep -i -E "genero_proximo|especie_candidata" /tmp/string_species_v12.txt
# Resultado: espécie mais próxima (ex: Asparagus officinalis) NÃO está disponível
# no STRING v12 → decisão de usar múltiplas espécies de referência (ortologia),
# seguindo o precedente de Carvalho et al. (2023).

# Obter os taxon IDs das 5 espécies de referência usadas no artigo:
grep -i -E "arabidopsis thaliana|musa acuminata|oryza sativa$|sorghum bicolor|zea mays" \
    /tmp/string_species_v12.txt
grep -i "oryza sativa" /tmp/string_species_v12.txt   # escolher a variedade "core" (Japonica)

# Taxon IDs confirmados:
#   Arabidopsis thaliana  = 3702
#   Sorghum bicolor       = 4558
#   Zea mays               = 4577
#   Musa acuminata         = 214687
#   Oryza sativa Japonica  = 39947

# =============================================================================
# ETAPA 3 — Baixar os proteomas das 5 espécies de referência (já em formato STRING)
# =============================================================================
mkdir -p "${ORTHO_DIR}/proteomas_referencia"
cd "${ORTHO_DIR}/proteomas_referencia"

for taxid in 3702 4558 4577 214687 39947; do
    curl -s -o "${taxid}.fa.gz" \
        "https://stringdb-downloads.org/download/protein.sequences.v12.0/${taxid}.protein.sequences.v12.0.fa.gz"
done
gunzip *.gz
ls -lh

# Conferir formato do header (esperado: ">taxid.proteinID")
for f in *.fa; do echo "== $f =="; grep ">" "$f" | head -n 2; done

# =============================================================================
# ETAPA 4 — Filtrar o proteoma da espécie de estudo para conter só os DEGs
# =============================================================================
# 4a) Extrair lista de IDs de DEGs (padj < 0.05 e |log2FC| > 1) da tabela DESeq2
awk -F'\t' 'NR>1 && $5!="NA" && $5<0.05 && ($2>1 || $2<-1) {print $6}' \
    "${DEG_TABLE}" > "${ORTHO_DIR}/lista_ids_deg.txt"
wc -l "${ORTHO_DIR}/lista_ids_deg.txt"

# 4b) Filtrar o .pep para manter só as sequências dos DEGs
awk 'BEGIN{
        while((getline id < "'"${ORTHO_DIR}"'/lista_ids_deg.txt")>0) degs[id]=1
     }
     /^>/{split($1,a,"."); a[1]=substr(a[1],2); keep=(a[1] in degs)}
     {if(keep) print}' \
    "${PEP_FILE}" > "${ORTHO_DIR}/proteomas_referencia/especie_estudo_degs.fa"

grep -c ">" "${ORTHO_DIR}/proteomas_referencia/especie_estudo_degs.fa"

# 4c) Validação cruzada do filtro (método alternativo, deve bater com 4b)
grep ">" "${PEP_FILE}" | sed -E 's/^>(\S+)\.p[0-9]+.*/\1/' | sort -u \
    > /tmp/pep_transcript_ids.txt
comm -12 <(sort "${ORTHO_DIR}/lista_ids_deg.txt") /tmp/pep_transcript_ids.txt | wc -l

# =============================================================================
# ETAPA 5 — Rodar OrthoFinder (espécie de estudo + 5 referências)
# =============================================================================
conda activate orthofinder_env

cd "${ORTHO_DIR}"
nohup orthofinder -f proteomas_referencia -t 16 -a 4 \
    > orthofinder.log 2>&1 &
disown
echo "PID: $!"

# Acompanhar progresso:
tail -n 30 "${ORTHO_DIR}/orthofinder.log"

# Resultado esperado ao final: resumo com nº de genes atribuídos a ortogrupos,
# em ${ORTHO_DIR}/proteomas_referencia/OrthoFinder/Results_<data>/

# =============================================================================
# ETAPA 6 — Cruzar ortogrupos: espécie de estudo (DEGs) x espécies de referência
# =============================================================================
pip install pandas --quiet   # se necessário

RESULTS_DIR="${ORTHO_DIR}/proteomas_referencia/OrthoFinder/Results_XXXXX"  # ajustar data

python3 << 'EOF'
import pandas as pd

RESULTS_DIR = "[diretório_rede_ortologia]/proteomas_referencia/OrthoFinder/Results_XXXXX"  # ajustar no script real
DEG_TABLE   = "[diretório_analise_deg]/DEGs_resultado_completo.tsv"
OUT_MAPPING = "[diretório_rede_ortologia]/tabela_mapeamento_ortologia.tsv"

# Direção dos DEGs (up/down), a partir da tabela completa do DESeq2
deg = pd.read_csv(DEG_TABLE, sep="\t")
deg = deg[(deg["padj"].notna()) & (deg["padj"] < 0.05) & (deg["log2FoldChange"].abs() > 1)]
deg["direction"] = deg["log2FoldChange"].apply(lambda x: "up" if x > 1 else "down")
direction_map = dict(zip(deg["gene"], deg["direction"]))

# Ortogrupos (colunas: Orthogroup, <taxid1..5>, especie_estudo_degs)
df = pd.read_csv(f"{RESULTS_DIR}/Orthogroups/Orthogroups.tsv", sep="\t").fillna("")
ref_cols = ["3702", "4558", "4577", "214687", "39947"]
study_col = "especie_estudo_degs"

rows = []
for _, row in df.iterrows():
    study_cell = row[study_col].strip()
    if not study_cell:
        continue
    ref_ids = []
    for col in ref_cols:
        if row[col].strip():
            ref_ids.extend(x.strip() for x in row[col].split(","))
    if not ref_ids:
        continue
    for ag in study_cell.split(","):
        ag = ag.strip()
        transcript_id = ag.rsplit(".p", 1)[0]
        direction = direction_map.get(transcript_id, "unknown")
        for rid in ref_ids:
            rows.append((transcript_id, direction, rid, row["Orthogroup"]))

out = pd.DataFrame(rows, columns=["transcrito_estudo", "direction", "string_ref_id", "orthogroup"])
out.to_csv(OUT_MAPPING, sep="\t", index=False)
print("Linhas na tabela final:", len(out))
EOF

# =============================================================================
# ETAPA 7 — Consultar a API do STRING para obter interações das espécies de referência
# =============================================================================
python3 << 'EOF'
import pandas as pd
import requests
import time

MAPPING_FILE = "[diretório_rede_ortologia]/tabela_mapeamento_ortologia.tsv"
OUT_RAW = "[diretório_rede_ortologia]/string_interacoes_brutas.tsv"

df = pd.read_csv(MAPPING_FILE, sep="\t")
df["taxid"] = df["string_ref_id"].str.split(".").str[0]

all_results = []
for taxid in df["taxid"].unique():
    ids = df[df["taxid"] == taxid]["string_ref_id"].unique()
    for i in range(0, len(ids), 1000):
        batch = ids[i:i+1000]
        r = requests.post(
            "https://string-db.org/api/tsv/network",
            data={
                "identifiers": "%0d".join(batch),
                "species": int(taxid),
                "required_score": 700,  # corte de confiança 0.70, igual ao artigo de referência
            },
        )
        if r.status_code == 200 and r.text.strip():
            lines = r.text.strip().split("\n")
            if len(lines) > 1:
                all_results.append(r.text)
        time.sleep(1)  # boa prática: não sobrecarregar a API do STRING

header = None
combined_lines = []
for text in all_results:
    lines = text.strip().split("\n")
    if header is None:
        header = lines[0]
        combined_lines.append(header)
    combined_lines.extend(lines[1:])

with open(OUT_RAW, "w") as f:
    f.write("\n".join(combined_lines))

print("Total de linhas de interação:", len(combined_lines) - 1)
EOF

# =============================================================================
# ETAPA 8 — Mapear as interações de volta para os transcritos da espécie de estudo
# =============================================================================
python3 << 'EOF'
import pandas as pd

MAPPING_FILE = "[diretório_rede_ortologia]/tabela_mapeamento_ortologia.tsv"
RAW_INTERACTIONS = "[diretório_rede_ortologia]/string_interacoes_brutas.tsv"
OUT_NETWORK = "[diretório_rede_ortologia]/rede_ppi_final.tsv"

mapping = pd.read_csv(MAPPING_FILE, sep="\t")
interactions = pd.read_csv(RAW_INTERACTIONS, sep="\t")

ref_to_study = mapping.groupby("string_ref_id").apply(
    lambda g: list(zip(g["transcrito_estudo"], g["direction"])), include_groups=False
).to_dict()

edges = []
for _, row in interactions.iterrows():
    a_list = ref_to_study.get(row["stringId_A"], [])
    b_list = ref_to_study.get(row["stringId_B"], [])
    for (ag_a, dir_a) in a_list:
        for (ag_b, dir_b) in b_list:
            if ag_a != ag_b:
                edges.append((ag_a, dir_a, ag_b, dir_b, row["score"]))

edf = pd.DataFrame(edges, columns=["no_A", "direction_A", "no_B", "direction_B", "score"])
edf["pair_key"] = edf.apply(lambda r: tuple(sorted([r["no_A"], r["no_B"]])), axis=1)
edf = edf.sort_values("score", ascending=False).drop_duplicates("pair_key")
edf = edf.drop(columns="pair_key")

edf.to_csv(OUT_NETWORK, sep="\t", index=False)
print("Total de arestas únicas na rede final:", len(edf))
print(edf[["direction_A", "direction_B"]].value_counts())
EOF

# =============================================================================
# FIM — Rede PPI final pronta em: ${ORTHO_DIR}/rede_ppi_final.tsv
# Colunas: no_A, direction_A, no_B, direction_B, score
# Próximos passos (não incluídos aqui): métricas topológicas (networkx / Network
# Analyzer no Cytoscape) e identificação de módulos (MCODE), conforme Carvalho
# et al. (2023).
# =============================================================================
