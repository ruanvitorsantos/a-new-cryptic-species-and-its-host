#!/usr/bin/env bash
# NCBI nr database, DIAMOND database (SLURM job) and a DIAMOND search of the assembly.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


# locate/download the BLAST nr volumes
find / -name update_blastdb.pl 2>/dev/null
[diretório_miniconda]/bin/update_blastdb.pl --decompress nr
# for BLAST the database flag is simply:  -db nr

# convert nr to a DIAMOND database (SLURM)
sbatch --job-name=nr_to_diamond \
  --partition=fila256 \
  --cpus-per-task=128 \
  --mem=900G \
  --time=7-00:00:00 \
  --wrap="export PATH=[diretório_miniconda]/bin:\$PATH && blastdbcmd -db [diretório_banco_nr]/nr -entry all -outfmt %f | diamond makedb --in - -d [diretório_banco_nr]/nr --threads 128"

# DIAMOND blastx of the assembly against nr (inside a SLURM job: threads = $SLURM_CPUS_PER_TASK)
diamond blastx \
  --query trinity_plantae.fasta \
  --db [diretório_banco_nr]/nr.dmnd \
  --out plantae_diamond_nr_outfmt6.tsv \
  --outfmt 6 \
  --max-target-seqs 5 \
  --evalue 1e-5 \
  --threads $SLURM_CPUS_PER_TASK

# generic local BLAST database
# makeblastdb -in dado.fna -dbtype nucl -out db_dado -parse_seqids
# blastn -query dado.fasta -db db_dado -task blastn -outfmt 6 -out resultado_dado.tsv -evalue 1e-10 -num_threads 8
