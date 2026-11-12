# Full analysis (6 diseased + 5 healthy): tximport at transcript level, low-count filter (filterByExpr), DESeq2, VST and a first PCA. Saves dds_completo.rds and vsd_completo.rds used by the next scripts.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

library(tximport)
library(DESeq2)
library(edgeR)
library(ggplot2)

dir_quant <- "[diretório_salmon_quants]"
dir_out <- "[diretório_deseq2_results]"
dir.create(dir_out, showWarnings = FALSE)

samples <- data.frame(
  sample = c("ASD2","ASD3","ASD5","ASD6","ASD7","ASD8",
             "ASS2","ASS5","ASS7","ASS8","ASS9"),
  condition = c(rep("Doente", 6), rep("Sadia", 5))
)

files <- setNames(
  file.path(dir_quant, c(
    "clean_ASD2_concatenado_S28_S17_S6_matched",
    "clean_ASD3_concatenado_S29_S18_S7_matched",
    "clean_ASD5_concatenado_S30_S19_S8_matched",
    "clean_ASD6_concatenado_S31_S20_S9_matched",
    "clean_ASD7_concatenado_S32_S21_S10_matched",
    "clean_ASD8_concatenado_S33_S22_S11_matched",
    "clean_ASS2_concatenado_S23_S12_S1_matched",
    "clean_ASS5_concatenado_S24_S13_S2_matched",
    "clean_ASS7_concatenado_S25_S14_S3_matched",
    "clean_ASS8_concatenado_S26_S15_S4_matched",
    "clean_ASS9_concatenado_S27_S16_S5_matched"
  ), "quant.sf"),
  samples$sample
)

txi <- tximport(files, type = "salmon", txOut = TRUE)

samples$condition <- factor(samples$condition, levels = c("Sadia", "Doente"))
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~condition)

keep <- filterByExpr(counts(dds), group = samples$condition)
dds <- dds[keep, ]
cat("Transcritos após filtro:", nrow(dds), "\n")

dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3) +
  scale_color_manual(values = c("Sadia" = "blue", "Doente" = "red")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

ggsave(file.path(dir_out, "PCA_inicial.png"), plot = p, width = 8, height = 6, dpi = 300)
saveRDS(dds, file.path(dir_out, "dds_completo.rds"))
saveRDS(vsd, file.path(dir_out, "vsd_completo.rds"))
cat("Concluído.\n")
