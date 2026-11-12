# Cross-check of the DESeq2 result with limma-voom (same filtered universe; exploratory).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("tximport", "edgeR", "limma"))

library(tximport)
library(edgeR)
library(limma)

dir_salmon <- "[diretório_salmon_quants]"

amostras <- c(
  "ASD2_concatenado_S28_S17_S6_matched", "ASD3_concatenado_S29_S18_S7_matched",
  "ASD5_concatenado_S30_S19_S8_matched", "ASD6_concatenado_S31_S20_S9_matched",
  "ASD7_concatenado_S32_S21_S10_matched", "ASD8_concatenado_S33_S22_S11_matched",
  "ASS2_concatenado_S23_S12_S1_matched", "ASS5_concatenado_S24_S13_S2_matched",
  "ASS7_concatenado_S25_S14_S3_matched", "ASS8_concatenado_S26_S15_S4_matched",
  "ASS9_concatenado_S27_S16_S5_matched"
)
condicao <- factor(c(rep("Doente", 6), rep("Sadia", 5)))

arquivos <- file.path(dir_salmon, amostras, "quant.sf")
names(arquivos) <- amostras
txi <- tximport(arquivos, type = "salmon", txOut = TRUE, dropInfReps = TRUE)

# filter the full txi first, so that the DGEList has the same universe as DESeq2
dge_temp <- DGEList(counts = txi$counts, group = condicao)
keep <- filterByExpr(dge_temp)
txi$counts    <- txi$counts[keep, ]
txi$abundance <- txi$abundance[keep, ]
txi$length    <- txi$length[keep, ]

dge <- DGEList(counts = txi$counts, group = condicao)
dge <- calcNormFactors(dge)          # TMM
cat("Transcritos após filtro:", nrow(dge), "\n")

design <- model.matrix(~ 0 + condicao)
colnames(design) <- levels(condicao)
v <- voom(dge, design, plot = TRUE)

fit <- lmFit(v, design)
cont.matrix <- makeContrasts(ASD_vs_ASS = Doente - Sadia, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont.matrix))

resLimma <- topTable(fit2, coef = "ASD_vs_ASS", number = Inf, sort.by = "P")
resUp   <- subset(resLimma, adj.P.Val < 0.05 & logFC >  2)    # logFC > 2 = fold change > 4
resDown <- subset(resLimma, adj.P.Val < 0.05 & logFC < -2)
cat("UP em Doente (ASD):", nrow(resUp), "\n")
cat("UP em Sadia  (ASS):", nrow(resDown), "\n")
