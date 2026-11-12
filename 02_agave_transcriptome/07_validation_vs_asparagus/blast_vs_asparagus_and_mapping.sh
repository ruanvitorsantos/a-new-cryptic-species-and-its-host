#!/usr/bin/env bash
# Similarity of the A. sisalana transcriptome to Asparagus officinalis (BLASTn) and read-mapping check of 100 random transcripts (bowtie2/samtools/IGV).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


# --- BLASTn against the A. officinalis genome ---
makeblastdb -in GCA_001876935.Asparagus_officinalis.fna -dbtype nucl -out db_asparagus -parse_seqids
blastn -query Agave_sisalana_transcriptome_awk.fasta -db db_asparagus -task blastn \
       -outfmt 6 -out resultado_agave_vs_aspargo.tsv -evalue 1e-10 -num_threads 8

# identity/mismatch summary of the alignments
python3 summarize_blast_identity.py resultado_agave_vs_aspargo.tsv

# --- mapping of reads to 100 random transcripts ---
conda install -c bioconda seqtk bowtie2 samtools -y
seqtk sample -s42 Agave_sisalana_transcriptome_awk.fasta 100 > 100_transcritos_aleatorios.fasta

bowtie2-build 100_transcritos_aleatorios.fasta bowtie2/indice_sisal_100
bowtie2 -x bowtie2/indice_sisal_100 -1 sua_amostra_R1.fastq.gz -2 sua_amostra_R2.fastq.gz \
        -S alinhamento_amostra_1.sam -p 8

samtools sort -@ 6 alinhamento_amostra_1.sam -o alinhamento_amostra_1_sorted.bam
samtools index alinhamento_amostra_1_sorted.bam           # .bai needed by IGV

# per-transcript depth and number of reads (input .bam/.bai + the 100-transcript FASTA are loaded in IGV)
samtools coverage alinhamento_amostra_1_sorted.bam -o estatisticas_cobertura.tsv
