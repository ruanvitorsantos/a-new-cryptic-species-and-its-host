# Venn diagram of the DEG sets of groups A and B (numbers entered from the DESeq2 results of each group).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
library(VennDiagram)

pdf("diagrama_venn.pdf", width = 7, height = 7)
grid.newpage()
draw.pairwise.venn(
  area1 = 43803,       # total group A
  area2 = 18861,       # total group B
  cross.area = 16260,  # intersection (A and B)
  category = c("Grupo A", "Grupo B"),
  fill = c("red", "blue"),
  alpha = c(0.5, 0.5),
  lty = "blank",
  cex = 1.3, cat.cex = 1.3,
  cat.pos = c(-20, 20), cat.dist = c(0.05, 0.05)
)
dev.off()
