# Directory placeholders

Every local path was replaced by a placeholder so that the scripts are the same for every user. Replace the placeholder with your own path.

| Placeholder | Meaning |
| --- | --- |
| `[diretório_projeto_sisal]` | Root directory of the *Agave sisalana* project |
| `[diretório_transcriptoma]` | Current transcriptome work directory: TransDecoder proteins (`protein/`), STRING files (`string_db/`, `string_db5/`), DIAMOND results (`ref5/`, `ref7/`) and the outputs of the network analysis |
| `[diretório_transcriptoma_base]` | Initial transcriptome work directory: assembly, rRNA databases, SortMeRNA results and BUSCO lineages |
| `[diretório_deseq2_results]` | Directory where the DESeq2 objects and tables are written (`dds_completo.rds`, `vsd_completo.rds`, `DEGs_todos_shrink.tsv`, ...) |
| `[diretório_tabelas_suplementares]` | Output directory of the supplementary tables |
| `[diretório_salmon]` | Salmon work directory (cleaned reads and indices) |
| `[diretório_salmon_quants]` | Directory with one sub-folder per sample containing `quant.sf` |
| `[diretório_interproscan]` | InterProScan installation directory (contains `interproscan.sh` and `resultados/`) |
| `[diretório_analise_deg]` | Directory with the DESeq2 result tables (used by the preliminary PPI script only) |
| `[diretório_rede_ortologia]` | Work directory of the preliminary orthology-based PPI script |
| `[diretório_banco_nr]` | Directory with the NCBI nr database (BLAST and DIAMOND formats) |
| `[diretório_banco_nt]` | Directory with the NCBI nt database |
| `[diretório_busco_aspergillus_odb12]` | BUSCO `aspergillus_odb12` lineage dataset |
| `[diretório_proteinas_busco]` | Output directory for the single-copy BUSCO proteins (one `.faa` per genome) |
| `[diretório_filogenia]` | Directory with the protein sets given to OrthoFinder |
| `[diretório_marcadores_cam]` | Directory with the marker (calmodulin) protein database |
| `[diretório_miniconda]` | Miniconda installation directory |

Conda/micromamba environments are referenced by **name** (for example `micromamba activate ncbi_datasets`), not by path.
