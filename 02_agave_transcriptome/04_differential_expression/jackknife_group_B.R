# Leave-one-out (jackknife) robustness of the DEGs, group B (healthy + diseased samples of group B).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

# Each DESeq2 run drops one sample; robust DEGs = intersection of all runs; Jaccard similarity heatmap between runs.

# 1. Carregar Bibliotecas Necessárias
library(DESeq2)
library(tximport)
library(tidyverse)
library(pheatmap)

# 2. Configuração de Diretórios
dir_base <- "[diretório_projeto_sisal]"
dir_quants <- file.path(dir_base, "salmon/salmon_quants")
dir_out <- file.path(dir_base, "salmon/jackknife_results")

# ==============================================================================
# 3. DEFINIÇÃO DAS AMOSTRAS DO GRUPO B (8 Amostras no Total)
# ==============================================================================
amostras_grupoB <- data.frame(
  sample = c(
    # 5 Controles (ASS)
    "ASS2_concatenado_S23_S12_S1_matched",
    "ASS5_concatenado_S24_S13_S2_matched",
    "ASS7_concatenado_S25_S14_S3_matched",
    "ASS8_concatenado_S26_S15_S4_matched",
    "ASS9_concatenado_S27_S16_S5_matched",
    
    # 3 Infectados do Grupo B (ASD 3, 7 e 8)
    "ASD3_concatenado_S29_S18_S7_matched",
    "ASD7_concatenado_S32_S21_S10_matched",
    "ASD8_concatenado_S33_S22_S11_matched"
  ),
  condition = factor(c(rep("Control", 5), rep("Infected", 3)))
)
amostras_grupoB$condition <- relevel(amostras_grupoB$condition, ref = "Control")

# 4. Mapeamento e Validação dos Arquivos
files_B <- file.path(dir_quants, amostras_grupoB$sample, "quant.sf")
names(files_B) <- amostras_grupoB$sample

if(!all(file.exists(files_B))) {
  cat("\nERRO CRÍTICO: Os seguintes arquivos não foram encontrados:\n")
  print(files_B[!file.exists(files_B)])
  stop("O script parou. Verifique os nomes no diretório.")
}

# 5. Funções Internas do Script
rodar_deseq2 <- function(txi_sub, amostras_sub) {
  dds <- DESeqDataSetFromTximport(txi = txi_sub, colData = amostras_sub, design = ~ condition)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds <- DESeq(dds, quiet = TRUE)
  
  res <- results(dds, contrast = c("condition", "Infected", "Control"))
  res_shrunk <- lfcShrink(dds, coef = "condition_Infected_vs_Control", type = "ashr", quiet = TRUE)
  
  df <- as.data.frame(res_shrunk) %>%
    rownames_to_column("transcript") %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 1)
  
  return(df)
}

calcular_jaccard <- function(set1, set2) { 
  length(intersect(set1, set2)) / length(union(set1, set2)) 
}

# 6. Importação de Dados
cat("\nImportando dados do Grupo B...\n")
txi_B <- tximport(files_B, type = "salmon", txOut = TRUE, countsFromAbundance = "lengthScaledTPM")

# 7. Análise Base do Grupo B
cat("\nRodando DESeq2 Base...\n")
base_degs_B <- rodar_deseq2(txi_B, amostras_grupoB)
cat("DEGs encontrados na base do Grupo B:", nrow(base_degs_B), "\n")

# 8. Pipeline Jackknife (Leave-One-Out)
jackknife_list_B <- list()
nomes_amostras_B <- amostras_grupoB$sample

cat("\nIniciando validação Jackknife para o Grupo B...\n")
for (amostra_excluida in nomes_amostras_B) {
  cat("  -> Removendo amostra:", amostra_excluida, "\n")
  
  amostras_sub_B <- amostras_grupoB[amostras_grupoB$sample != amostra_excluida, ]
  txi_sub_B <- list(
    abundance = txi_B$abundance[, amostras_sub_B$sample],
    counts = txi_B$counts[, amostras_sub_B$sample],
    length = txi_B$length[, amostras_sub_B$sample],
    countsFromAbundance = "lengthScaledTPM"
  )
  
  degs_parciais_B <- rodar_deseq2(txi_sub_B, amostras_sub_B)
  jackknife_list_B[[paste0("Sem_", amostra_excluida)]] <- degs_parciais_B$transcript
}

# 9. Interseção e Exportação de DEGs Robustos
todos_conjuntos_B <- c(list(Base = base_degs_B$transcript), jackknife_list_B)
degs_robustos_B <- Reduce(intersect, todos_conjuntos_B)

df_robustos_B <- base_degs_B %>% filter(transcript %in% degs_robustos_B)
write_delim(df_robustos_B, file.path(dir_out, "DEGs_GrupoB_Robustos.tsv"), delim = "\t")

cat("\n======================================================\n")
cat("RESULTADOS FINAIS - GRUPO B\n")
cat("Total de DEGs Robustos (Verdadeiros):", length(degs_robustos_B), "\n")
cat("Taxa de Falsos Positivos identificada:", nrow(base_degs_B) - length(degs_robustos_B), "transcritos.\n")
cat("======================================================\n")

# 10. Matriz de Jaccard e Heatmap
jaccard_matrix_B <- matrix(NA, length(todos_conjuntos_B), length(todos_conjuntos_B))
rownames(jaccard_matrix_B) <- colnames(jaccard_matrix_B) <- names(todos_conjuntos_B)

for(i in 1:length(todos_conjuntos_B)) {
  for(j in 1:length(todos_conjuntos_B)) {
    jaccard_matrix_B[i,j] <- calcular_jaccard(todos_conjuntos_B[[i]], todos_conjuntos_B[[j]])
  }
}

rownames(jaccard_matrix_B) <- gsub("_concatenado.*", "", rownames(jaccard_matrix_B))
colnames(jaccard_matrix_B) <- gsub("_concatenado.*", "", colnames(jaccard_matrix_B))

arquivo_pdf <- file.path(dir_out, "Heatmap_Outlier_Jackknife_GrupoB.pdf")
pdf(arquivo_pdf, width = 9, height = 7)
pheatmap(jaccard_matrix_B, 
         display_numbers = TRUE, 
         color = colorRampPalette(c("white", "gold", "firebrick3"))(50),
         main = "Similaridade (Jaccard) - Grupo B",
         fontsize_number = 10, angle_col = "45")
invisible(dev.off())

cat("\nHeatmap salvo com sucesso em:", arquivo_pdf, "\n")
