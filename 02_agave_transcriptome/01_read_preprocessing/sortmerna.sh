#!/usr/bin/env bash
# Removal of ribosomal RNA with SortMeRNA (Section 5.2.1).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

# The default database (smr_v4.3_default_db.fasta) covers eukaryotic 18S/28S, bacterial and archaeal 16S/23S, and 5S/5.8S (Rfam).

conda install -c bioconda sortmerna=4.3.4 -y
git clone https://github.com/biocore/sortmerna.git      # source of the rRNA databases

sortmerna \
  --ref [diretório_transcriptoma]/rRNA_databases/smr_v4.3_default_db.fasta \
  --reads [diretório_transcriptoma]/sortmerna_results/Agave_sisalana_transcriptome.fasta \
  --workdir [diretório_transcriptoma]/sortmerna_results/sortmerna_out \
  --aligned reads_rRNA \
  --other Agave_sisalana_transcriptoma_sem_rRNA \
  --fastx \
  --threads 4
