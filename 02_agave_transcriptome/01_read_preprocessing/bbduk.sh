#!/usr/bin/env bash
# Adapter removal, quality trimming (Phred 30) and k-mer enrichment with BBDuk/BBTools (Section 5.2.1).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


conda install -c bioconda bbmap -y

# (1) adapter removal + quality trimming in one pass
bbduk.sh in1=amostra_R1.fastq.gz in2=amostra_R2.fastq.gz \
         out1=amostra_R1_processada.fastq.gz out2=amostra_R2_processada.fastq.gz \
         ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo \
         qtrim=rl trimq=30

# (the same in two separate steps)
# bbduk.sh in1=amostra_R1.fastq.gz in2=amostra_R2.fastq.gz out1=amostra_R1_no_adapter.fastq.gz out2=amostra_R2_no_adapter.fastq.gz ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo
# bbduk.sh in1=amostra_R1_no_adapter.fastq.gz in2=amostra_R2_no_adapter.fastq.gz out1=amostra_R1_limpa.fastq.gz out2=amostra_R2_limpa.fastq.gz qtrim=rl trimq=30

# (2) k-mer matching against a reference FASTA: reads that match go to outm, the rest to out
#     (used with the reference genomes of related Agave species to enrich plant reads)
bbduk.sh in1=reads_R1.fastq in2=reads_R2.fastq ref=sua_referencia.fasta \
         outm=reads_alvo_encontrado.fastq out=reads_limpas.fastq k=31 stats=estatisticas_bbduk.txt

# (3) summary table from the *_stats.txt files (output directory: bbduk_filter_agave/)
echo -e "Amostra\tTotal_Reads\tReads_Planta\t(%)\tReads_Alvo_Unmatched\t(%)" > bbduk_filter_agave/resumo_agave_bbduk.tsv
for stats in bbduk_filter_agave/*_stats.txt; do
    nome=$(basename "$stats" _stats.txt)
    total=$(grep -w "Input:" "$stats" | awk '{print $2}')
    planta=$(grep "Total Removed:" "$stats" | awk '{print $3}')
    pct_planta=$(grep "Total Removed:" "$stats" | awk '{print $5}' | tr -d '()')
    alvo=$(grep -w "Result:" "$stats" | awk '{print $2}')
    pct_alvo=$(grep -w "Result:" "$stats" | awk '{print $4}' | tr -d '()')
    echo -e "${nome}\t${total}\t${planta}\t${pct_planta}\t${alvo}\t${pct_alvo}" >> bbduk_filter_agave/resumo_agave_bbduk.tsv
done
column -t -s $'\t' bbduk_filter_agave/resumo_agave_bbduk.tsv
