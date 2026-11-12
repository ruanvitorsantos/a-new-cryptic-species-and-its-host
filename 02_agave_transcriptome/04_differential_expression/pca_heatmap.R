# Exploratory analyses (Section 5.2.4): VST, PCA and heatmap of the 500 most variable transcripts.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

# Requires the `dds` object created by deseq2.R (load("deseq2_workspace.RData")).

library(DESeq2)
library(pheatmap)
library(matrixStats)
library(viridisLite)

# VST (blind = FALSE)
vsd <- vst(dds, blind = FALSE)

# PCA: plotPCA uses the 500 transcripts with the highest variance by default
pdf("PCA.pdf", width = 6, height = 6)
plotPCA(vsd, intgroup = "condition")
dev.off()

# Heatmap of the top 500 most variable transcripts (Z-score per transcript)
rv     <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[1:500]
mat    <- assay(vsd)[select, ]

mat_scaled <- t(scale(t(mat)))
mat_scaled[is.na(mat_scaled)]       <- 0      # remove NA/Inf produced by zero-variance rows
mat_scaled[is.infinite(mat_scaled)] <- 0

annotation_col <- data.frame(Condition = colData(vsd)$condition)
rownames(annotation_col) <- colnames(mat_scaled)

pheatmap(mat_scaled,
         cluster_rows = TRUE, cluster_cols = FALSE,
         color = viridisLite::viridis(100),
         show_rownames = FALSE, show_colnames = TRUE,
         annotation_col = annotation_col,
         fontsize_col = 8, border_color = NA,
         filename = "heatmap_top500.pdf", width = 8, height = 10)
