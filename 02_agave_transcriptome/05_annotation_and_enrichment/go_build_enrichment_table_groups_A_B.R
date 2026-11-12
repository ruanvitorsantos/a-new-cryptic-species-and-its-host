# Assembles Enriquecimento_GO_Completo.tsv (groups A and B, Up/Down, BP/MF/CC), the input of go_plot_nonredundant_groups_A_B.R and of Table S7.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

# [CONFIRMAR] Reconstructed from the table layout and from go_enrichment_groups_A_B.R / go_add_ontology_class.R / go_restore_ids_and_round.R; the original
# script that wrote this table was not kept. Terms are keyed by GO ID here (term names are attached with TERM2NAME), which gives the same statistics.

library(tidyverse)
library(clusterProfiler)
library(GO.db)
library(AnnotationDbi)

dir_base <- "[diretório_projeto_sisal]"
universe_ids <- readLines(file.path(dir_base, "salmon/universe_ids.txt"))

# DEGs of each group (padj < 0.05 and |log2FC| > 1)
extrair_degs <- function(arquivo, up = TRUE) {
  df <- read_delim(file.path(dir_base, "salmon", arquivo), delim = "\t", show_col_types = FALSE)
  if (up) {
    df %>% filter(padj < 0.05 & log2FoldChange > 1) %>% pull(1)
  } else {
    df %>% filter(padj < 0.05 & log2FoldChange < -1) %>% pull(1)
  }
}

# GO dictionary from the InterProScan TSV (column 14), transcript ids without the TransDecoder .pN suffix
interpro_raw <- read_delim(file.path(dir_base, "interproscan-5.76-107.0/resultados/ASS_anotacao.tsv"),
                           delim = "\t", col_names = paste0("V", 1:15), show_col_types = FALSE)
go_df <- interpro_raw %>%
  mutate(transcript_id = gsub("\\.p[0-9]+", "", V1)) %>%
  dplyr::select(V14, transcript_id) %>%
  filter(V14 != "-") %>%
  separate_rows(V14, sep = "\\|") %>%
  mutate(go_id = str_extract(V14, "GO:[0-9]+")) %>%
  filter(!is.na(go_id)) %>%
  distinct(go_id, transcript_id)

info <- suppressMessages(AnnotationDbi::select(GO.db, keys = unique(go_df$go_id),
                                               columns = c("TERM", "ONTOLOGY"), keytype = "GOID")) %>%
  distinct(GOID, .keep_all = TRUE)
term2name <- info %>% transmute(go_id = GOID, name = str_to_title(ifelse(is.na(TERM), GOID, TERM)))

# Over-representation test (hypergeometric, Benjamini-Hochberg), one call per group and direction
run_enrich <- function(genes, grupo, regulacao) {
  res <- enricher(gene = genes, universe = universe_ids, TERM2GENE = go_df, TERM2NAME = term2name,
                  pAdjustMethod = "BH", pvalueCutoff = 0.05)
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
  as.data.frame(res) %>% mutate(Grupo = grupo, Regulacao = regulacao, GO_ID = ID)
}

classe_pt <- c(BP = "Processo Biológico", MF = "Função Molecular", CC = "Componente Celular")

tabela <- bind_rows(
  run_enrich(extrair_degs("DEGs_DEG_A_shrink.tsv", TRUE),  "Grupo A", "Up"),
  run_enrich(extrair_degs("DEGs_DEG_A_shrink.tsv", FALSE), "Grupo A", "Down"),
  run_enrich(extrair_degs("DEGs_DEG_B_shrink.tsv", TRUE),  "Grupo B", "Up"),
  run_enrich(extrair_degs("DEGs_DEG_B_shrink.tsv", FALSE), "Grupo B", "Down")
) %>%
  left_join(info %>% transmute(GO_ID = GOID, Classe = unname(classe_pt[ONTOLOGY])), by = "GO_ID") %>%
  mutate(FoldEnrichment = round(as.numeric(FoldEnrichment), 2),
         pvalue = signif(as.numeric(pvalue), 3),
         p.adjust = signif(as.numeric(p.adjust), 3),
         qvalue = signif(as.numeric(qvalue), 3)) %>%
  dplyr::select(Grupo, Regulacao, Classe, GO_ID, Description, FoldEnrichment, GeneRatio, BgRatio,
                pvalue, p.adjust, qvalue, Count, geneID) %>%
  arrange(Grupo, Regulacao, Classe, p.adjust)

write_delim(tabela, file.path(dir_base, "salmon/Enriquecimento_GO_Completo.tsv"), delim = "\t")
cat("Table written:", nrow(tabela), "enriched terms\n")
