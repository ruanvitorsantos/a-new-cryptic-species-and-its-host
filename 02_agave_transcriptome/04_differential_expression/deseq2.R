# Differential expression at the transcript level with DESeq2 (Section 5.2.3), contrast D (Doente/diseased) vs H (Sadia/healthy).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

# Reference level: Sadia (healthy). Samples: 6 ASD (diseased) + 5 ASS (healthy).

library(tximport)
library(DESeq2)
library(edgeR)      # filterByExpr
library(apeglm)

dir_salmon <- "[diretório_salmon_quants]"

amostras <- c("ASD2_concatenado_S28_S17_S6_matched", "ASD3_concatenado_S29_S18_S7_matched",
              "ASD5_concatenado_S30_S19_S8_matched", "ASD6_concatenado_S31_S20_S9_matched",
              "ASD7_concatenado_S32_S21_S10_matched", "ASD8_concatenado_S33_S22_S11_matched",
              "ASS2_concatenado_S23_S12_S1_matched", "ASS5_concatenado_S24_S13_S2_matched",
              "ASS7_concatenado_S25_S14_S3_matched", "ASS8_concatenado_S26_S15_S4_matched",
              "ASS9_concatenado_S27_S16_S5_matched")

condicao <- factor(c(rep("Doente", 6), rep("Sadia", 5)), levels = c("Sadia", "Doente"))   # Sadia = reference
coldata  <- data.frame(row.names = amostras, condition = condicao)

arquivos <- file.path(dir_salmon, amostras, "quant.sf")
names(arquivos) <- amostras

# txOut = TRUE keeps the analysis at the transcript level (required for Trinity assemblies)
txi <- tximport(arquivos, type = "salmon", txOut = TRUE, dropInfReps = TRUE)
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)

# [CONFIRMAR] low-count filter taken from Methods 5.2.3; it is not in the original notes
keep <- filterByExpr(counts(dds), group = coldata$condition)
dds  <- dds[keep, ]

dds <- DESeq(dds)

# [CONFIRMAR] unshrunken results and apeglm-shrunken log2FC, as in Methods 5.2.3 (the notes only show the unshrunken step)
res        <- results(dds, contrast = c("condition", "Doente", "Sadia"))
res_shrink <- lfcShrink(dds, coef = "condition_Doente_vs_Sadia", type = "apeglm")

resOrdered <- res_shrink[order(res_shrink$pvalue), ]

# DEGs: padj < 0.05 and |log2FC| > 1
resUp   <- subset(resOrdered, padj < 0.05 & log2FoldChange >  1)   # up in diseased (ASD)
resDown <- subset(resOrdered, padj < 0.05 & log2FoldChange < -1)   # [CONFIRMAR] the notes had < -2 here, while the comment and Methods say |log2FC| > 1

cat("Up:", nrow(resUp), " Down:", nrow(resDown), " Total:", nrow(resUp) + nrow(resDown), "\n")

write.table(as.data.frame(res_shrink), "DEGs_shrink.tsv", sep = "\t", quote = FALSE)
save.image("deseq2_workspace.RData")
