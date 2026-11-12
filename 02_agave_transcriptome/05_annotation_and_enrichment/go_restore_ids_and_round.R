# Utility: restore GO identifiers and round numeric columns of the complete enrichment table.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


library(tidyverse)
library(GO.db)
library(AnnotationDbi)

cat("\n--- RESTAURANDO CÓDIGOS GO E AJUSTANDO NÚMEROS ---\n")

# 1. Carregar a planilha atual
caminho_planilha <- "[diretório_salmon]/Enriquecimento_GO_Completo.tsv"
tabela <- read_delim(caminho_planilha, delim = "\t", show_col_types = FALSE)

# 2. Recriar o dicionário exato a partir do InterPro
interpro_raw <- read_delim("[diretório_interproscan]/resultados/ASS_anotacao.tsv", delim = "\t", col_names = paste0("V", 1:15), show_col_types = FALSE)

go_df_clean <- interpro_raw %>%
  dplyr::select(V14) %>%
  filter(V14 != "-") %>%
  separate_rows(V14, sep = "\\|") %>%
  mutate(go_id = str_extract(V14, "GO:[0-9]+")) %>%
  filter(!is.na(go_id)) %>%
  distinct(go_id)

termos <- suppressMessages(AnnotationDbi::select(GO.db, keys = go_df_clean$go_id, columns = "TERM", keytype = "GOID"))

# Criar o de/para com a mesma formatação
de_para <- data.frame(
  go_id = go_df_clean$go_id,
  term = termos$TERM[match(go_df_clean$go_id, termos$GOID)],
  stringsAsFactors = FALSE
) %>%
  mutate(
    term = ifelse(is.na(term), go_id, term),
    term = str_to_title(term)
  ) %>%
  distinct(term, .keep_all = TRUE)

# 3. Cruzar, corrigir a ID e arredondar para o Excel
tabela_limpa <- tabela %>%
  left_join(de_para, by = c("Description" = "term")) %>%
  mutate(
    ID = ifelse(!is.na(go_id), go_id, ID),
    FoldEnrichment = round(as.numeric(FoldEnrichment), 2),
    zScore = round(as.numeric(zScore), 2),
    pvalue = signif(as.numeric(pvalue), 3),
    p.adjust = signif(as.numeric(p.adjust), 3),
    qvalue = signif(as.numeric(qvalue), 3)
  ) %>%
  dplyr::select(-go_id)

# 4. Sobrescrever arquivo
write_delim(tabela_limpa, caminho_planilha, delim = "\t")

cat("Sucesso! Coluna ID atualizada para os códigos GO e colunas matemáticas arredondadas.\n")
