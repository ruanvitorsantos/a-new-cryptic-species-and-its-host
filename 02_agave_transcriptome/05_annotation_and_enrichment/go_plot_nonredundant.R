# Non-redundant GO dotplots: terms sharing >= 80% of genes (overlap coefficient, cutoff = 0.8) are collapsed, keeping the most significant;
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

# the top 25 terms per direction (Up/Down) are plotted for each group and GO class (BP/CC/MF).
# Input: Enriquecimento_GO_Completo.tsv (all enrichment results, with columns Grupo, Classe, Regulacao, ID/GO_ID, Description, Count, geneID, p.adjust).

library(tidyverse)

# 1. Carregar e padronizar dados
dir_base <- "[diretório_projeto_sisal]"
tabela_plot <- read_delim(file.path(dir_base, "salmon/Enriquecimento_GO_Completo.tsv"), delim = "\t", show_col_types = FALSE) %>%
  mutate(
    Codigo_GO = if("GO_ID" %in% names(.)) GO_ID else if(any(grepl("GO:", ID))) ID else NA,
    Rotulo_Completo = paste0(Codigo_GO, " ", str_to_title(Description))
  ) %>%
  filter(!is.na(Classe))

# 2. Filtro de Redundância Otimizado
limpar_redundancia <- function(df, cutoff = 0.8) {
  if(nrow(df) <= 1) return(df)
  df <- arrange(df, p.adjust)
  
  mantidos <- c()
  genes_list <- list()
  
  for (i in seq_len(nrow(df))) {
    g_atuais <- unlist(strsplit(df$geneID[i], "/"))
    
    # Compara o overlap com todos os genes já guardados de forma vetorizada
    if (length(genes_list) > 0) {
      overlaps <- sapply(genes_list, function(g_comp) length(intersect(g_atuais, g_comp)) / min(length(g_atuais), length(g_comp)))
      if (any(overlaps >= cutoff)) next
    }
    
    mantidos <- c(mantidos, df$Rotulo_Completo[i])
    genes_list[[length(genes_list) + 1]] <- g_atuais
  }
  
  filter(df, Rotulo_Completo %in% mantidos)
}

# 3. Função de Plotagem Compacta
gerar_plot <- function(grupo_alvo, classe_alvo) {
  df <- tabela_plot %>% filter(Grupo == grupo_alvo, Classe == classe_alvo)
  if (nrow(df) == 0) return(invisible())
  
  df_top <- df %>%
    group_by(Regulacao) %>%
    group_modify(~ limpar_redundancia(.x, 0.8)) %>%
    slice_max(Count, n = 25, with_ties = FALSE) %>%
    arrange(Regulacao, Count) %>%
    ungroup() %>%
    mutate(Rotulo_Completo = factor(Rotulo_Completo, levels = unique(Rotulo_Completo)))
  
  p <- ggplot(df_top, aes(x = Count, y = Rotulo_Completo)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    facet_grid(Regulacao ~ ., scales = "free_y", space = "free_y") +
    labs(title = paste("Enriquecimento GO:", grupo_alvo, "-", classe_alvo), x = "Gene Count", y = NULL) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(face = "bold")
    ) +
    scale_color_gradientn(colours = c("#FFD700", "#008000", "#1B4D3E"), breaks = c(0, 0.025, 0.05), guide = guide_colorbar(reverse = TRUE))
  
  sigla <- c("Processo Biológico"="BP", "Componente Celular"="CC", "Função Molecular"="MF")[classe_alvo]
  arquivo <- file.path(dir_base, sprintf("salmon/Figura_%s_%s_Limpo.pdf", gsub(" ", "", grupo_alvo), sigla))
  
  ggsave(arquivo, plot = p, width = 15, height = 10)
}

# 4. Execução em Lote Automatizada
expand_grid(
  grupo = c("Grupo A", "Grupo B"), 
  classe = c("Processo Biológico", "Componente Celular", "Função Molecular")
) %>%
  pwalk(~ gerar_plot(.x, .y))
