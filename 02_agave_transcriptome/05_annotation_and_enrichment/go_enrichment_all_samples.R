# GO over-representation analysis (clusterProfiler::enricher, BH, FDR < 0.05) of up- and down-regulated DEGs of the full analysis; universe = all tested transcripts; GO annotation from InterProScan.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

library(tidyverse)
library(clusterProfiler)
library(GO.db)
library(AnnotationDbi)

dir_out <- "[diretório_deseq2_results]"
arq_interpro <- "[diretório_interproscan]/resultados/ASS_novo_anotacao.tsv"

res_df <- read.table(file.path(dir_out, "DEGs_todos_shrink.tsv"), header=TRUE, sep="\t")
universe_ids <- res_df$gene
deg_up <- res_df %>% dplyr::filter(padj < 0.05 & log2FoldChange > 1) %>% dplyr::pull(gene)
deg_down <- res_df %>% dplyr::filter(padj < 0.05 & log2FoldChange < -1) %>% dplyr::pull(gene)

interpro_raw <- read_delim(arq_interpro, delim="\t",
                           col_names=paste0("V", 1:15),
                           show_col_types=FALSE)

go_df_clean <- interpro_raw %>%
  mutate(transcript_id = gsub("\\.p[0-9]+", "", V1)) %>%
  dplyr::select(V14, transcript_id) %>%
  dplyr::filter(V14 != "-") %>%
  separate_rows(V14, sep="\\|") %>%
  mutate(go_id = str_extract(V14, "GO:[0-9]+")) %>%
  dplyr::filter(!is.na(go_id)) %>%
  distinct(go_id, transcript_id)

go_df_id <- data.frame(
  gs_name = go_df_clean$go_id,
  gene = go_df_clean$transcript_id,
  stringsAsFactors=FALSE
)

termos <- suppressMessages(AnnotationDbi::select(GO.db,
          keys=unique(go_df_clean$go_id),
          columns=c("TERM","ONTOLOGY"), keytype="GOID"))

term2name <- data.frame(
  go_id = termos$GOID,
  term = str_to_title(termos$TERM),
  stringsAsFactors=FALSE
)

run_enrich <- function(genes) {
  enricher(gene=genes, universe=universe_ids,
           TERM2GENE=go_df_id,
           TERM2NAME=term2name,
           pAdjustMethod="BH", pvalueCutoff=0.05)
}

go_up <- run_enrich(deg_up)
go_down <- run_enrich(deg_down)
cat("Termos UP:", nrow(as.data.frame(go_up)), "\n")
cat("Termos DOWN:", nrow(as.data.frame(go_down)), "\n")

adicionar_classe <- function(df, regulacao) {
  ont_map <- setNames(termos$ONTOLOGY, termos$GOID)
  df$Classe <- dplyr::recode(ont_map[df$ID],
                      "BP"="Processo Biológico",
                      "MF"="Função Molecular",
                      "CC"="Componente Celular")
  df$Regulacao <- regulacao
  return(df)
}

df_up <- adicionar_classe(as.data.frame(go_up), "Up")
df_down <- adicionar_classe(as.data.frame(go_down), "Down")
tabela_completa <- bind_rows(df_up, df_down) %>%
  dplyr::filter(Classe %in% c("Processo Biológico", "Função Molecular")) %>%
  arrange(Classe, Regulacao, p.adjust)

write.table(tabela_completa,
            file.path(dir_out, "Enriquecimento_BP_MF_completo.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

saveRDS(go_up, file.path(dir_out, "go_up.rds"))
saveRDS(go_down, file.path(dir_out, "go_down.rds"))

cat("BP e MF - Up:", nrow(dplyr::filter(tabela_completa, Regulacao=="Up")), "\n")
cat("BP e MF - Down:", nrow(dplyr::filter(tabela_completa, Regulacao=="Down")), "\n")
cat("Concluído.\n")
