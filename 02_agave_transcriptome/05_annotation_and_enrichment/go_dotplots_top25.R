# GO dotplots (BP and MF): the 25 most significant terms (smallest adjusted p-value) for up- and down-regulated transcripts.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

library(tidyverse)
library(ggplot2)

dir_out <- "[diretório_deseq2_results]"
tabela <- read.table(file.path(dir_out, "Enriquecimento_BP_MF_completo.tsv"),
                     header=TRUE, sep="\t")

tema_artigo <- theme_bw(base_size=13) +
  theme(
    panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold", size=15, hjust=0.5),
    axis.title=element_text(face="bold"),
    axis.text.y=element_text(color="black", size=10),
    strip.background=element_rect(fill="gray95", color="black"),
    strip.text=element_text(face="bold", size=12)
  )

gerar_plot <- function(classe_alvo, arquivo_saida) {
  df <- tabela %>%
    dplyr::filter(Classe == classe_alvo) %>%
    dplyr::group_by(Regulacao) %>%
    dplyr::slice_min(p.adjust, n=25, with_ties=FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Count = as.numeric(sapply(strsplit(GeneRatio, "/"), `[`, 1)),
      Rotulo = paste0(ID, " ", str_wrap(Description, width=45))
    ) %>%
    dplyr::arrange(Regulacao, Count) %>%
    dplyr::mutate(
      order_id = row_number(),
      label = paste0(Regulacao, "_", order_id, "_", Rotulo),
      label = factor(label, levels=unique(label))
    )

  if (nrow(df) == 0) {
    cat("Sem dados para:", classe_alvo, "\n")
    return(invisible())
  }

  p <- ggplot(df, aes(x=Count, y=label)) +
    geom_point(aes(size=Count, color=p.adjust)) +
    facet_grid(Regulacao ~ ., scales="free_y", space="free_y") +
    scale_y_discrete(labels=function(x) gsub("^[^_]+_[0-9]+_", "", x)) +
    labs(title=paste("Enriquecimento GO:", classe_alvo),
         x="Gene Count", y=NULL) +
    tema_artigo +
    scale_color_gradientn(colours=c("#FFD700","#008000","#1B4D3E"),
                          breaks=c(0, 0.025, 0.05),
                          guide=guide_colorbar(reverse=TRUE)) +
    scale_size_continuous(range=c(4,10))

  ggsave(arquivo_saida, plot=p, width=17, height=11)
  cat("Salvo:", arquivo_saida, "\n")
}

gerar_plot("Processo Biológico", file.path(dir_out, "Enriquecimento_BP_v4.pdf"))
gerar_plot("Função Molecular", file.path(dir_out, "Enriquecimento_MF_v4.pdf"))
cat("Concluído.\n")
