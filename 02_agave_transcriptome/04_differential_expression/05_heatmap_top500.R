# Heatmap of the 500 most variable transcripts (VST, Z-score per transcript, viridis).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

library(DESeq2)
library(pheatmap)
library(matrixStats)
library(viridisLite)
library(dplyr)

dir_out <- "[diretório_deseq2_results]"
vsd <- readRDS(file.path(dir_out, "vsd_completo.rds"))

rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[1:500]
mat <- assay(vsd)[select, ]
mat_scaled <- t(scale(t(mat)))
mat_scaled[is.na(mat_scaled)] <- 0
mat_scaled[is.infinite(mat_scaled)] <- 0

annotation_col <- data.frame(Condition = colData(vsd)$condition)
rownames(annotation_col) <- colnames(mat_scaled)

png(file.path(dir_out, "Heatmap_top500.png"), width = 10, height = 12,
    units = "in", res = 300)
pheatmap(mat_scaled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = viridisLite::viridis(100),
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = annotation_col,
         fontsize_col = 8,
         border_color = NA)
dev.off()
cat("Concluído.\n")
