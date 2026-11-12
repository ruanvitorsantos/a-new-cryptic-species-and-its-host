# MA plot of the apeglm-shrunken results.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

library(ggplot2)
library(dplyr)

dir_out <- "[diretório_deseq2_results]"
res_df <- read.table(file.path(dir_out, "DEGs_todos_shrink.tsv"),
                     header=TRUE, sep="\t")

res_df <- res_df %>%
  dplyr::filter(!is.na(padj) & !is.na(log2FoldChange) & baseMean > 0) %>%
  dplyr::mutate(
    A = log2(baseMean),
    M = log2FoldChange,
    grupo = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    ),
    grupo = factor(grupo, levels=c("NS","Down","Up"))
  ) %>%
  dplyr::arrange(grupo)

p <- ggplot(res_df, aes(x=A, y=M, color=grupo)) +
  geom_point(size=0.8, alpha=0.5) +
  scale_color_manual(values=c("Up"="red","Down"="blue","NS"="gray80"),
                     name="Expressão") +
  geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=0.8) +
  geom_hline(yintercept=c(-1,1), linetype="dotted", color="black", alpha=0.7) +
  coord_cartesian(ylim=c(-8, 8)) +
  labs(title="MA Plot - ASD vs ASS",
       x="A (log2 Mean Expression)",
       y="M (log2 Fold Change)") +
  theme_bw(base_size=13) +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        panel.grid.minor=element_blank())

ggsave(file.path(dir_out, "MA_plot_v2.pdf"), plot=p, width=10, height=7)
cat("Concluído.\n")
