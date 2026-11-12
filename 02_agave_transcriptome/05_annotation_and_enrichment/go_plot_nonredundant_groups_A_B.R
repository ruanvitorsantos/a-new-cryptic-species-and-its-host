# Non-redundant GO dotplots for groups A and B: terms sharing >= 80% of genes are collapsed (overlap coefficient), keeping the most significant; top 25 per direction.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

library(tidyverse)
library(ggplot2)

cat("\n--- INICIANDO PIPELINE COM FILTRO DE SIMILARIDADE (0.8) E N=25 ---\n")

dir_base <- "[diretório_projeto_sisal]"
tabela <- read_delim(file.path(dir_base, "salmon/Enriquecimento_GO_Completo.tsv"), delim = "\t", show_col_types = FALSE)

tabela_plot <- tabela %>%
  mutate(
    Codigo_GO = if("GO_ID" %in% names(.)) GO_ID else if(any(grepl("GO:", ID))) ID else NA,
    Rotulo_Completo = paste0(Codigo_GO, " ", str_to_title(Description))
  ) %>%
  filter(!is.na(Classe))

remover_redundancia <- function(df_sub, cutoff = 0.8) {
  if(nrow(df_sub) <= 1) return(df_sub)
  df_sub <- df_sub %>% arrange(p.adjust)
  termos_mantidos <- c()
  lista_genes_mantidos <- list()
  for (i in 1:nrow(df_sub)) {
    genes_atuais <- unlist(strsplit(df_sub$geneID[i], "/"))
    redundante <- FALSE
    if (length(lista_genes_mantidos) > 0) {
      for (j in 1:length(lista_genes_mantidos)) {
        genes_comparacao <- lista_genes_mantidos[[j]]
        intersecao <- length(intersect(genes_atuais, genes_comparacao))
        overlap <- intersecao / min(length(genes_atuais), length(genes_comparacao))
        if (overlap >= cutoff) {
          redundante <- TRUE
          break
        }
      }
    }
    if (!redundante) {
      termos_mantidos <- c(termos_mantidos, df_sub$Rotulo_Completo[i])
      lista_genes_mantidos[[length(lista_genes_mantidos) + 1]] <- genes_atuais
    }
  }
  return(df_sub %>% filter(Rotulo_Completo %in% termos_mantidos))
}

tema <- theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black", size = 10),
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(color = "gray95"),
    strip.background = element_rect(fill = "gray95", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing.y = unit(0.3, "lines")
  )

plotar_padronizado_limpo <- function(df, grupo_alvo, classe_alvo, arquivo_saida) {
  df_filtro <- df %>% filter(Grupo == grupo_alvo, Classe == classe_alvo)
  if (nrow(df_filtro) == 0) {
    cat("Sem enriquecimento significativo para:", grupo_alvo, "-", classe_alvo, "\n")
    return(NULL)
  }
  df_limpo <- df_filtro %>%
    group_by(Regulacao) %>%
    group_modify(~ remover_redundancia(.x, cutoff = 0.8)) %>%
    ungroup()
  df_top <- df_limpo %>%
    group_by(Regulacao) %>%
    slice_max(Count, n = 25, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(Regulacao, Count) %>%
    mutate(Rotulo_Completo = factor(Rotulo_Completo, levels = unique(Rotulo_Completo)))
  p <- ggplot(df_top, aes(x = Count, y = Rotulo_Completo)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    facet_grid(Regulacao ~ ., scales = "free_y") +
    labs(title = paste("Enriquecimento GO:", grupo_alvo, "-", classe_alvo), x = "Gene Count", y = NULL) +
    tema +
    scale_color_gradientn(colours = c("#FFD700", "#008000", "#1B4D3E"), breaks = c(0, 0.025, 0.05), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title = "p-adjust", reverse = TRUE)) +
    scale_size_continuous(range = c(3, 8), guide = guide_legend(title = "Gene Count", override.aes = list(color = "black")))
  ggsave(arquivo_saida, plot = p, width = 15, height = 10, units = "in", dpi = 300)
  cat("Ficheiro salvo:", arquivo_saida, "\n")
}

cat("\nA processar figuras do Grupo A...\n")
plotar_padronizado_limpo(tabela_plot, "Grupo A", "Processo Biológico", file.path(dir_base, "salmon/Figura_GrupoA_BP_Limpo.pdf"))
plotar_padronizado_limpo(tabela_plot, "Grupo A", "Componente Celular", file.path(dir_base, "salmon/Figura_GrupoA_CC_Limpo.pdf"))
plotar_padronizado_limpo(tabela_plot, "Grupo A", "Função Molecular", file.path(dir_base, "salmon/Figura_GrupoA_MF_Limpo.pdf"))

cat("\nA processar figuras do Grupo B...\n")
plotar_padronizado_limpo(tabela_plot, "Grupo B", "Processo Biológico", file.path(dir_base, "salmon/Figura_GrupoB_BP_Limpo.pdf"))
plotar_padronizado_limpo(tabela_plot, "Grupo B", "Componente Celular", file.path(dir_base, "salmon/Figura_GrupoB_CC_Limpo.pdf"))
plotar_padronizado_limpo(tabela_plot, "Grupo B", "Função Molecular", file.path(dir_base, "salmon/Figura_GrupoB_MF_Limpo.pdf"))
