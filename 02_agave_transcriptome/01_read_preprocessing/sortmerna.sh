#!/usr/bin/env bash
# Removal of ribosomal RNA with SortMeRNA (Section 5.2.1).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

# The default database (smr_v4.3_default_db.fasta) covers eukaryotic 18S/28S, bacterial and archaeal 16S/23S, and 5S/5.8S (Rfam).

# SortMeRNA v8.0.1 was used (installed from the official release installer).
# [CONFIRMAR] the command below was recorded with an earlier SortMeRNA version; check the flags and the input (reads vs. assembly) for v8.0.1.
git clone https://github.com/biocore/sortmerna.git      # source of the rRNA databases

sortmerna \
  --ref [diretório_transcriptoma]/rRNA_databases/smr_v4.3_default_db.fasta \
  --reads [diretório_transcriptoma]/sortmerna_results/Agave_sisalana_transcriptome.fasta \
  --workdir [diretório_transcriptoma]/sortmerna_results/sortmerna_out \
  --aligned reads_rRNA \
  --other Agave_sisalana_transcriptoma_sem_rRNA \
  --fastx \
  --threads 4
