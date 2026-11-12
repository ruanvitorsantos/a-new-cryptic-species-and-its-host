# Volcano plot (apeglm-shrunken log2FC). `dds_A` is the DESeq2 object of the analysis to be plotted (e.g. `dds` from deseq2.R).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


library(DESeq2)
library(ggplot2)

res_A <- results(dds_A, contrast = c("condition", "Doente", "Sadia"))
res_A_shrink <- lfcShrink(dds_A, coef = "condition_Doente_vs_Sadia", type = "apeglm")

res_A_df <- na.omit(as.data.frame(res_A_shrink))

res_A_df$group <- "NS"
res_A_df$group[res_A_df$padj < 0.05 & res_A_df$log2FoldChange >  1] <- "ASD"   # up in diseased
res_A_df$group[res_A_df$padj < 0.05 & res_A_df$log2FoldChange < -1] <- "ASS"   # up in healthy

pdf("Volcano_A.pdf", width = 6, height = 6)
ggplot(res_A_df, aes(x = log2FoldChange, y = -log10(padj), color = group)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("ASD" = "red", "ASS" = "blue", "NS" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano Plot - DEG-A", x = "log2 Fold Change", y = "-log10 adjusted p-value") +
  theme_minimal()
dev.off()
