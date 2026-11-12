# apeglm-shrunken log2FC and DEG selection (padj < 0.05, |log2FC| > 1) from dds_completo.rds; writes DEGs_todos_shrink.tsv (basis of Table S6).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

library(tximport)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(dplyr)

dir_out <- "[diretório_deseq2_results]"
dds <- readRDS(file.path(dir_out, "dds_completo.rds"))

# Extrair resultados com shrinkage
res <- lfcShrink(dds, coef = "condition_Doente_vs_Sadia", type = "apeglm")
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Filtrar DEGs significativos
deg <- res_df %>% filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
cat("Total DEGs:", nrow(deg), "\n")
cat("Up-regulated:", nrow(filter(deg, log2FoldChange > 1)), "\n")
cat("Down-regulated:", nrow(filter(deg, log2FoldChange < -1)), "\n")

# Salvar
write.table(res_df, file.path(dir_out, "DEGs_todos_shrink.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat("Concluído.\n")
