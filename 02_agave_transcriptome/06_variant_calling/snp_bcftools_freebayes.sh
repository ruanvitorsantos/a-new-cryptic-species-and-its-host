#!/usr/bin/env bash
# SNP calling on a single transcript of interest (diploid bcftools calling and 5n FreeBayes calling).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

# Requires a sorted/indexed BAM of the reads aligned to the assembly (see ../07_validation_vs_asparagus for the bowtie2/samtools step).

conda create -n snp_analysis -c bioconda -c conda-forge bcftools samtools gsl -y
conda activate snp_analysis
conda install -c bioconda freebayes -y

# extract the transcript of interest from the original BAM
samtools view -b alinhamento_amostra_1_sorted.bam TRINITY_DN591067_c1_g2_i3 > alvo.bam
samtools index alvo.bam
samtools faidx trinity_plantae.fasta      # .fai index of the reference

# --- bcftools (diploid) ---
bcftools mpileup -d 1000 -f trinity_plantae.fasta alvo.bam | \
  bcftools call -mv -Ov -o snp_diploide.vcf
bcftools filter -i 'QUAL > 30 && DP > 50' snp_diploide.vcf -o snp_final_bcftools.vcf
grep -v "^#" snp_diploide.vcf | wc -l          # number of variants

# --- FreeBayes (ploidy 5) ---
freebayes -f trinity_plantae.fasta -p 5 alvo.bam > freebayes_raw.vcf
bcftools filter -i 'QUAL > 30 && INFO/DP > 50' freebayes_raw.vcf -o snp_final_freebayes.vcf
