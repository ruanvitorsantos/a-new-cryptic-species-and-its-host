# Stratified DESeq2 analyses: group A (ASD2, ASD5, ASD6) and group B (ASD3, ASD7, ASD8), each against the five healthy samples.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

library(tximport)
library(DESeq2)
library(edgeR)
library(dplyr)

dir_quant <- "[diretório_salmon_quants]"
dir_out <- "[diretório_deseq2_results]"

todos_arquivos <- c(
  ASD2 = file.path(dir_quant, "clean_ASD2_concatenado_S28_S17_S6_matched"),
  ASD3 = file.path(dir_quant, "clean_ASD3_concatenado_S29_S18_S7_matched"),
  ASD5 = file.path(dir_quant, "clean_ASD5_concatenado_S30_S19_S8_matched"),
  ASD6 = file.path(dir_quant, "clean_ASD6_concatenado_S31_S20_S9_matched"),
  ASD7 = file.path(dir_quant, "clean_ASD7_concatenado_S32_S21_S10_matched"),
  ASD8 = file.path(dir_quant, "clean_ASD8_concatenado_S33_S22_S11_matched"),
  ASS2 = file.path(dir_quant, "clean_ASS2_concatenado_S23_S12_S1_matched"),
  ASS5 = file.path(dir_quant, "clean_ASS5_concatenado_S24_S13_S2_matched"),
  ASS7 = file.path(dir_quant, "clean_ASS7_concatenado_S25_S14_S3_matched"),
  ASS8 = file.path(dir_quant, "clean_ASS8_concatenado_S26_S15_S4_matched"),
  ASS9 = file.path(dir_quant, "clean_ASS9_concatenado_S27_S16_S5_matched")
)

rodar_grupo <- function(amostras, nome_grupo) {
  files <- paste0(todos_arquivos[amostras], "/quant.sf")
  names(files) <- amostras
  samples <- data.frame(
    sample = amostras,
    condition = factor(ifelse(grepl("ASD", amostras), "Doente", "Sadia"), levels=c("Sadia","Doente")),
    row.names = amostras
  )
  txi <- tximport(files, type="salmon", txOut=TRUE)
  dds <- DESeqDataSetFromTximport(txi, colData=samples, design=~condition)
  keep <- filterByExpr(counts(dds), group=samples$condition)
  dds <- dds[keep,]
  dds <- DESeq(dds)
  res <- lfcShrink(dds, coef="condition_Doente_vs_Sadia", type="apeglm")
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  deg <- res_df %>% filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
  cat(nome_grupo, "- Total DEGs:", nrow(deg), "\n")
  cat(nome_grupo, "- Up:", nrow(filter(deg, log2FoldChange > 1)), "\n")
  cat(nome_grupo, "- Down:", nrow(filter(deg, log2FoldChange < -1)), "\n")
  write.table(res_df, file.path(dir_out, paste0("DEGs_", nome_grupo, "_shrink.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  saveRDS(dds, file.path(dir_out, paste0("dds_", nome_grupo, ".rds")))
}

rodar_grupo(c("ASD2","ASD5","ASD6","ASS2","ASS5","ASS7","ASS8","ASS9"), "GrupoA")
rodar_grupo(c("ASD3","ASD7","ASD8","ASS2","ASS5","ASS7","ASS8","ASS9"), "GrupoB")
cat("Concluido.\n")
