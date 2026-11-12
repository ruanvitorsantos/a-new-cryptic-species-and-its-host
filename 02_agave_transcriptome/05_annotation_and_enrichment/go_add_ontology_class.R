# Adds the GO class (BP/MF/CC) to the enrichment results and writes Enriquecimento_BP_MF_completo.tsv.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

library(tidyverse)
library(GO.db)
library(AnnotationDbi)

dir_out <- "[diretório_deseq2_results]"

go_up <- readRDS(file.path(dir_out, "go_up.rds"))
go_down <- readRDS(file.path(dir_out, "go_down.rds"))

# Função para adicionar classe GO (BP, MF, CC)
adicionar_classe <- function(df, regulacao) {
  ids <- df$ID
  ids_limpos <- str_extract(ids, "GO:[0-9]+")
  ids_validos <- ids_limpos[!is.na(ids_limpos)]
  
  if (length(ids_validos) > 0) {
    ontologia <- suppressMessages(AnnotationDbi::select(GO.db,
                 keys=ids_validos, columns="ONTOLOGY", keytype="GOID"))
    
    classe_map <- setNames(ontologia$ONTOLOGY, ontologia$GOID)
    df$Classe <- ifelse(is.na(ids_limpos), NA, classe_map[ids_limpos])
    df$Classe <- recode(df$Classe,
                        "BP" = "Processo Biológico",
                        "MF" = "Função Molecular",
                        "CC" = "Componente Celular")
  } else {
    df$Classe <- NA
  }
  df$Regulacao <- regulacao
  return(df)
}

df_up <- adicionar_classe(as.data.frame(go_up), "Up")
df_down <- adicionar_classe(as.data.frame(go_down), "Down")
tabela_completa <- bind_rows(df_up, df_down)

# Filtrar apenas BP e MF
tabela_filtrada <- tabela_completa %>%
  filter(Classe %in% c("Processo Biológico", "Função Molecular")) %>%
  arrange(Classe, Regulacao, p.adjust)

write.table(tabela_filtrada,
            file.path(dir_out, "Enriquecimento_BP_MF_completo.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

cat("BP e MF - Up:", nrow(filter(tabela_filtrada, Regulacao=="Up")), "\n")
cat("BP e MF - Down:", nrow(filter(tabela_filtrada, Regulacao=="Down")), "\n")
cat("Planilha salva.\n")
