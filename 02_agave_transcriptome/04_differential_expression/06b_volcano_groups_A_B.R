# Volcano plots of the stratified analyses (groups A and B).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

library(ggplot2)
library(dplyr)

dir_out <- "[diretório_deseq2_results]"

gerar_volcano <- function(arquivo, titulo, nome_saida) {
  res_df <- read.table(arquivo, header=TRUE, sep="\t")
  res_df$group <- "NS"
  res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "ASD"
  res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "ASS"
  p <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=group)) +
    geom_point(size=1.5, alpha=0.7) +
    scale_color_manual(values=c("ASD"="red","ASS"="blue","NS"="gray")) +
    geom_vline(xintercept=c(-1,1), linetype="dashed") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    labs(title=titulo, x="log2 Fold Change", y="-log10 adjusted p-value") +
    theme_minimal()
  ggsave(file.path(dir_out, nome_saida), plot=p, width=8, height=6, dpi=300)
  cat(titulo, "- concluído.\n")
}

gerar_volcano(file.path(dir_out, "DEGs_GrupoA_shrink.tsv"), "Volcano - Grupo A", "Volcano_GrupoA.png")
gerar_volcano(file.path(dir_out, "DEGs_GrupoB_shrink.tsv"), "Volcano - Grupo B", "Volcano_GrupoB.png")
