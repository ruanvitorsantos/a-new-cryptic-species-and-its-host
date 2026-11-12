# Utility: look up the name and ontology class of a GO identifier with GO.db.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


# ==============================================================================
# SCRIPT 2: FERRAMENTA DE CONSULTA INDIVIDUAL DE ONTOLOGIA
# ==============================================================================

library(GO.db)
library(AnnotationDbi)

consultar_go <- function(codigo_go) {
  resultado <- suppressMessages(AnnotationDbi::select(GO.db, keys = codigo_go, columns = c("TERM", "ONTOLOGY"), keytype = "GOID"))
  
  classe <- switch(resultado$ONTOLOGY[1],
                   "BP" = "Processo Biológico",
                   "MF" = "Função Molecular",
                   "CC" = "Componente Celular",
                   "Desconhecido")
                   
  cat("\nCódigo:   ", resultado$GOID[1], 
      "\nNome:     ", resultado$TERM[1], 
      "\nClasse:   ", classe, "\n\n")
}

# Altere o código abaixo para consultar qualquer termo
consultar_go("GO:0005730")
