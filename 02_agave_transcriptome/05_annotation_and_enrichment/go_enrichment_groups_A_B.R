# GO over-representation analysis for the stratified analyses (groups A and B). Inputs: DEGs_DEG_A_shrink.tsv / DEGs_DEG_B_shrink.tsv, universe_ids.txt and the InterProScan TSV.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


# 1. Carregamento de Pacotes
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(GO.db)
library(AnnotationDbi)

# 2. Definição de Diretórios e Arquivos
dir_base <- "[diretório_projeto_sisal]"
arq_universe <- file.path(dir_base, "salmon/universe_ids.txt")
arq_interpro <- file.path(dir_base, "interproscan-5.76-107.0/resultados/ASS_anotacao.tsv")

# 3. Carregamento do Universo e DEGs (Grupos A e B)
universe_ids <- readLines(arq_universe)

# Função para extrair listas de DEGs (Filtro: padj < 0.05 e |log2FC| > 1)
extrair_degs <- function(caminho_arquivo, up = TRUE) {
  df <- read_delim(caminho_arquivo, delim = "\t", show_col_types = FALSE)
  if (up) {
    df %>% filter(padj < 0.05 & log2FoldChange > 1) %>% pull(1)
  } else {
    df %>% filter(padj < 0.05 & log2FoldChange < -1) %>% pull(1)
  }
}

deg_A_up <- extrair_degs(file.path(dir_base, "salmon/DEGs_DEG_A_shrink.tsv"), up = TRUE)
deg_A_down <- extrair_degs(file.path(dir_base, "salmon/DEGs_DEG_A_shrink.tsv"), up = FALSE)
deg_B_up <- extrair_degs(file.path(dir_base, "salmon/DEGs_DEG_B_shrink.tsv"), up = TRUE)
deg_B_down <- extrair_degs(file.path(dir_base, "salmon/DEGs_DEG_B_shrink.tsv"), up = FALSE)

# 4. Preparação e Limpeza do Dicionário InterProScan
interpro_raw <- read_delim(arq_interpro, delim = "\t", col_names = paste0("V", 1:15), show_col_types = FALSE)

go_df_clean <- interpro_raw %>%
  mutate(transcript_id = gsub("\\.p[0-9]+", "", V1)) %>%
  dplyr::select(V14, transcript_id) %>%
  filter(V14 != "-") %>%
  separate_rows(V14, sep = "\\|") %>%
  mutate(go_id = str_extract(V14, "GO:[0-9]+")) %>%
  filter(!is.na(go_id)) %>%
  distinct(go_id, transcript_id)

# 5. Tradução de Termos GO (Base R para evitar conflitos de pacotes)
termos <- suppressMessages(AnnotationDbi::select(GO.db, keys = unique(go_df_clean$go_id), columns = "TERM", keytype = "GOID"))

go_df_final <- data.frame(
  term = termos$TERM[match(go_df_clean$go_id, termos$GOID)],
  gene = go_df_clean$transcript_id,
  stringsAsFactors = FALSE
)

go_df_final$term <- ifelse(is.na(go_df_final$term), go_df_clean$go_id, go_df_final$term)
go_df_final$term <- str_to_title(go_df_final$term)

# 6. Testes de Enriquecimento (Teste Hipergeométrico)
run_enrich <- function(genes) {
  enricher(gene = genes, universe = universe_ids, TERM2GENE = go_df_final, pAdjustMethod = "BH", pvalueCutoff = 0.05)
}

go_A_up <- run_enrich(deg_A_up)
go_A_down <- run_enrich(deg_A_down)
go_B_up <- run_enrich(deg_B_up)
go_B_down <- run_enrich(deg_B_down)

# 7. Configuração Estética Global (ggplot2)
tema_artigo <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10, margin = margin(r = 5)),
    plot.title = element_text(face = "bold", hjust = 0.5),
    aspect.ratio = 1.3
  )

# 8. Função de Plotagem e Exportação
plotar_e_salvar <- function(go_obj, titulo, nome_arquivo) {
  if (!is.null(go_obj) && nrow(as.data.frame(go_obj)) > 0) {
    p <- dotplot(go_obj, showCategory = 20, label_format = 200) + 
      ggtitle(titulo) + 
      tema_artigo + 
      scale_color_viridis_c(option = "cividis", guide = guide_colorbar(reverse = TRUE))
    
    ggsave(file.path(dir_base, "salmon", nome_arquivo), plot = p, width = 13, height = 8)
  }
}

# 9. Geração dos Gráficos Finais
plotar_e_salvar(go_A_up, "GO - Ativados (Up) - Grupo A", "Figura_GO_A_Up.pdf")
plotar_e_salvar(go_A_down, "GO - Reprimidos (Down) - Grupo A", "Figura_GO_A_Down.pdf")
plotar_e_salvar(go_B_up, "GO - Ativados (Up) - Grupo B", "Figura_GO_B_Up.pdf")
plotar_e_salvar(go_B_down, "GO - Reprimidos (Down) - Grupo B", "Figura_GO_B_Down.pdf")
