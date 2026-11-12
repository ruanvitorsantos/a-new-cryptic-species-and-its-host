# Leave-one-out (jackknife) robustness of the DEGs, group A (5 healthy + 3 diseased).
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
dir.create(dir_out, showWarnings = FALSE)

# 3. DEFINIÇÃO DAS AMOSTRAS DO GRUPO A (8 Amostras no Total)
amostras_grupoA <- data.frame(
  sample = c(
    # 5 Controles (ASS)
    "ASS2_concatenado_S23_S12_S1_matched",
    "ASS5_concatenado_S24_S13_S2_matched",
    "ASS7_concatenado_S25_S14_S3_matched",
    "ASS8_concatenado_S26_S15_S4_matched",
    "ASS9_concatenado_S27_S16_S5_matched",
    
    # 3 Infectados do Grupo A (ASD 2, 5 e 6)
    "ASD2_concatenado_S28_S17_S6_matched",
    "ASD5_concatenado_S30_S19_S8_matched",
    "ASD6_concatenado_S31_S20_S9_matched"
  ),
  condition = factor(c(rep("Control", 5), rep("Infected", 3)))
)
amostras_grupoA$condition <- relevel(amostras_grupoA$condition, ref = "Control")

# 4. Mapeamento e Validação dos Arquivos
files_A <- file.path(dir_quants, amostras_grupoA$sample, "quant.sf")
names(files_A) <- amostras_grupoA$sample

if(!all(file.exists(files_A))) {
  cat("\nERRO CRÍTICO: Os seguintes arquivos não foram encontrados:\n")
  print(files_A[!file.exists(files_A)])
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
cat("\nImportando dados do Grupo A...\n")
txi_A <- tximport(files_A, type = "salmon", txOut = TRUE, countsFromAbundance = "lengthScaledTPM")

# 7. Análise Base do Grupo A
cat("\nRodando DESeq2 Base...\n")
base_degs_A <- rodar_deseq2(txi_A, amostras_grupoA)
cat("DEGs encontrados na base do Grupo A:", nrow(base_degs_A), "\n")

# 8. Pipeline Jackknife (Leave-One-Out)
jackknife_list_A <- list()
nomes_amostras_A <- amostras_grupoA$sample

cat("\nIniciando validação Jackknife para o Grupo A...\n")
for (amostra_excluida in nomes_amostras_A) {
  cat("  -> Removendo amostra:", amostra_excluida, "\n")
  
  amostras_sub_A <- amostras_grupoA[amostras_grupoA$sample != amostra_excluida, ]
  txi_sub_A <- list(
    abundance = txi_A$abundance[, amostras_sub_A$sample],
    counts = txi_A$counts[, amostras_sub_A$sample],
    length = txi_A$length[, amostras_sub_A$sample],
    countsFromAbundance = "lengthScaledTPM"
  )
  
  degs_parciais_A <- rodar_deseq2(txi_sub_A, amostras_sub_A)
  jackknife_list_A[[paste0("Sem_", amostra_excluida)]] <- degs_parciais_A$transcript
}

# 9. Interseção e Exportação de DEGs Robustos
todos_conjuntos_A <- c(list(Base = base_degs_A$transcript), jackknife_list_A)
degs_robustos_A <- Reduce(intersect, todos_conjuntos_A)

df_robustos_A <- base_degs_A %>% filter(transcript %in% degs_robustos_A)
write_delim(df_robustos_A, file.path(dir_out, "DEGs_GrupoA_Robustos_Final.tsv"), delim = "\t")
cat("RESULTADOS FINAIS - GRUPO A (Atualizado)\n")
cat("Total de DEGs Robustos (Verdadeiros):", length(degs_robustos_A), "\n")
cat("Taxa de Falsos Positivos identificada:", nrow(base_degs_A) - length(degs_robustos_A), "transcritos.\n")

# 10. Matriz de Jaccard e Heatmap
jaccard_matrix_A <- matrix(NA, length(todos_conjuntos_A), length(todos_conjuntos_A))
rownames(jaccard_matrix_A) <- colnames(jaccard_matrix_A) <- names(todos_conjuntos_A)

for(i in 1:length(todos_conjuntos_A)) {
  for(j in 1:length(todos_conjuntos_A)) {
    jaccard_matrix_A[i,j] <- calcular_jaccard(todos_conjuntos_A[[i]], todos_conjuntos_A[[j]])
  }
}

rownames(jaccard_matrix_A) <- gsub("_concatenado.*", "", rownames(jaccard_matrix_A))
colnames(jaccard_matrix_A) <- gsub("_concatenado.*", "", colnames(jaccard_matrix_A))

arquivo_pdf_A <- file.path(dir_out, "Heatmap_Outlier_Jackknife_GrupoA_Final.pdf")
pdf(arquivo_pdf_A, width = 9, height = 7)
pheatmap(jaccard_matrix_A, 
         display_numbers = TRUE, 
         color = colorRampPalette(c("white", "gold", "firebrick3"))(50),
         main = "Similaridade (Jaccard) - Grupo A",
         fontsize_number = 10, angle_col = "45")
invisible(dev.off())

cat("\nHeatmap salvo com sucesso em:", arquivo_pdf_A, "\n")
