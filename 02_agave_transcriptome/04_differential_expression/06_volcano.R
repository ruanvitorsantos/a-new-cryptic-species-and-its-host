# Volcano plot of the apeglm-shrunken results (diseased vs healthy).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

library(ggplot2)
library(dplyr)

dir_out <- "[diretório_deseq2_results]"
res_df <- read.table(file.path(dir_out, "DEGs_todos_shrink.tsv"),
                     header=TRUE, sep="\t")

res_df$group <- "NS"
res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up-regulated"
res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down-regulated"
res_df$group <- factor(res_df$group, levels=c("NS","Down-regulated","Up-regulated"))
res_df <- res_df %>% dplyr::arrange(group)

p <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=group)) +
  geom_point(size=1.5, alpha=0.7) +
  scale_color_manual(values=c("Up-regulated"="red",
                               "Down-regulated"="blue",
                               "NS"="gray80"),
                     name="Expressão") +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  labs(title="Volcano Plot - Agave sisalana infectada vs sadia",
       x="log2 Fold Change",
       y="-log10 adjusted p-value") +
  theme_minimal(base_size=13) +
  theme(plot.title=element_text(face="bold", hjust=0.5))

ggsave(file.path(dir_out, "Volcano_v2.pdf"), plot=p, width=9, height=7)
cat("Concluído.\n")
