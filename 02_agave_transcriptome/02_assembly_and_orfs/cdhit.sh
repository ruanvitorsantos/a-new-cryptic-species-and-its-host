#!/usr/bin/env bash
# Redundancy reduction / overlap between assemblies with CD-HIT.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).


# clustering of the assembly at 95% identity
cd-hit-est -i trinity_plantae.fasta -o similaridade_plantae_ch-hit -c 0.95 -n 10 -d 0 -M 0 -T 20

# comparison of the assembly of this study against an external A. sisalana reference transcriptome (90% identity)
cd-hit-est-2d -i Transcriptome_Codings_Agave_sisalana.fasta \
              -i2 Agave_sisalana_transcriptome.fasta \
              -o overlap_cdhit_90 \
              -c 0.9 -n 7 -d 0 -M 0 -T 20
