# Leave-one-out (jackknife) robustness of the DEGs, all 11 samples.
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
# 3. DEFINIÇÃO DE TODAS AS AMOSTRAS (11 Amostras no Total)
# ==============================================================================
amostras_todos <- data.frame(
  sample = c(
    # 5 Controles (ASS)
    "ASS2_concatenado_S23_S12_S1_matched",
    "ASS5_concatenado_S24_S13_S2_matched",
    "ASS7_concatenado_S25_S14_S3_matched",
    "ASS8_concatenado_S26_S15_S4_matched",
    "ASS9_concatenado_S27_S16_S5_matched",
    
    # 6 Infectados (ASD Grupo A + Grupo B)
    "ASD2_concatenado_S28_S17_S6_matched",
    "ASD5_concatenado_S30_S19_S8_matched",
    "ASD6_concatenado_S31_S20_S9_matched",
    "ASD3_concatenado_S29_S18_S7_matched",
    "ASD7_concatenado_S32_S21_S10_matched",
    "ASD8_concatenado_S33_S22_S11_matched"
  ),
  condition = factor(c(rep("Control", 5), rep("Infected", 6)))
)
amostras_todos$condition <- relevel(amostras_todos$condition, ref = "Control")

# 4. Mapeamento e Validação dos Arquivos
files_todos <- file.path(dir_quants, amostras_todos$sample, "quant.sf")
names(files_todos) <- amostras_todos$sample

if(!all(file.exists(files_todos))) {
  cat("\nERRO CRÍTICO: Arquivos ausentes:\n")
  print(files_todos[!file.exists(files_todos)])
  stop("O script parou. Verifique os diretórios.")
}

# 5. Função DESeq2
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

# 6. Importação
cat("\nImportando dados de todas as amostras...\n")
txi_todos <- tximport(files_todos, type = "salmon", txOut = TRUE, countsFromAbundance = "lengthScaledTPM")

# 7. Análise Base Unificada
cat("\nRodando DESeq2 Base Unificado...\n")
base_degs_todos <- rodar_deseq2(txi_todos, amostras_todos)
cat("DEGs encontrados na base unificada:", nrow(base_degs_todos), "\n")

# 8. Pipeline Jackknife
jackknife_list_todos <- list()
nomes_amostras_todos <- amostras_todos$sample

cat("\nIniciando validação Jackknife para todos os dados...\n")
for (amostra_excluida in nomes_amostras_todos) {
  cat("  -> Removendo amostra:", amostra_excluida, "\n")
  
  amostras_sub <- amostras_todos[amostras_todos$sample != amostra_excluida, ]
  txi_sub <- list(
    abundance = txi_todos$abundance[, amostras_sub$sample],
    counts = txi_todos$counts[, amostras_sub$sample],
    length = txi_todos$length[, amostras_sub$sample],
    countsFromAbundance = "lengthScaledTPM"
  )
  
  degs_parciais <- rodar_deseq2(txi_sub, amostras_sub)
  jackknife_list_todos[[paste0("Sem_", amostra_excluida)]] <- degs_parciais$transcript
}

# 9. Interseção e Exportação
todos_conjuntos <- c(list(Base = base_degs_todos$transcript), jackknife_list_todos)
degs_robustos_todos <- Reduce(intersect, todos_conjuntos)

df_robustos_todos <- base_degs_todos %>% filter(transcript %in% degs_robustos_todos)
write_delim(df_robustos_todos, file.path(dir_out, "DEGs_Todos_Robustos.tsv"), delim = "\t")

cat("\n======================================================\n")
cat("RESULTADOS FINAIS - ANÁLISE UNIFICADA\n")
cat("Total de DEGs Robustos (Verdadeiros):", length(degs_robustos_todos), "\n")
cat("Taxa de Falsos Positivos identificada:", nrow(base_degs_todos) - length(degs_robustos_todos), "transcritos.\n")
cat("======================================================\n")

# 10. Heatmap de Jaccard
jaccard_matrix_todos <- matrix(NA, length(todos_conjuntos), length(todos_conjuntos))
rownames(jaccard_matrix_todos) <- colnames(jaccard_matrix_todos) <- names(todos_conjuntos)

for(i in 1:length(todos_conjuntos)) {
  for(j in 1:length(todos_conjuntos)) {
    jaccard_matrix_todos[i,j] <- calcular_jaccard(todos_conjuntos[[i]], todos_conjuntos[[j]])
  }
}

rownames(jaccard_matrix_todos) <- gsub("_concatenado.*", "", rownames(jaccard_matrix_todos))
colnames(jaccard_matrix_todos) <- gsub("_concatenado.*", "", colnames(jaccard_matrix_todos))

arquivo_pdf <- file.path(dir_out, "Heatmap_Outlier_Jackknife_Todos.pdf")
pdf(arquivo_pdf, width = 11, height = 9)
pheatmap(jaccard_matrix_todos, 
         display_numbers = TRUE, 
         color = colorRampPalette(c("white", "gold", "firebrick3"))(50),
         main = "Similaridade (Jaccard) - Todos os Grupos",
         fontsize_number = 9, angle_col = "45")
invisible(dev.off())

cat("\nHeatmap salvo com sucesso em:", arquivo_pdf, "\n")
