# PCA of the VST data with the three sample groups outlined.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

library(DESeq2)
library(ggplot2)

dir_out <- "[diretório_deseq2_results]"
vsd <- readRDS(file.path(dir_out, "vsd_completo.rds"))

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_data$group <- ifelse(pca_data$name %in% c("ASD2","ASD5","ASD6"), "Group 1",
                  ifelse(pca_data$name %in% c("ASD3","ASD7","ASD8"), "Group 2",
                  "Control group"))

p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_vline(xintercept = 30, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  stat_ellipse(data = subset(pca_data, group == "Group 1"),
               aes(x = PC1, y = PC2), color = "black", linewidth = 0.8) +
  stat_ellipse(data = subset(pca_data, group == "Group 2"),
               aes(x = PC1, y = PC2), color = "black", linewidth = 0.8) +
  stat_ellipse(data = subset(pca_data, group == "Control group"),
               aes(x = PC1, y = PC2), color = "black", linewidth = 0.8) +
  geom_point(aes(color = condition), size = 4) +
  geom_text(aes(label = name, color = condition), vjust = -1, size = 3.5) +
  annotate("text", x = -80, y = 60, label = "Group 1", size = 5) +
  annotate("text", x = -50, y = -70, label = "Group 2", size = 5) +
  annotate("text", x = 75, y = 20, label = "Control group", size = 5) +
  scale_color_manual(values = c("Doente" = "#E69B2E", "Sadia" = "#4472C4"),
                     labels = c("Doente" = "Disease", "Sadia" = "Health"),
                     name = "Condition") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% de variance")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top left",
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"))

ggsave(file.path(dir_out, "PCA_novo_estetica.png"), plot = p,
       width = 10, height = 8, dpi = 300)
cat("Concluído.\n")
