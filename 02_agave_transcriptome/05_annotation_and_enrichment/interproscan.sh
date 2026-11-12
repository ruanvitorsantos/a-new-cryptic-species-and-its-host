#!/usr/bin/env bash
# Functional annotation of the TransDecoder proteins with InterProScan v5.76 (Section 5.2.5).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

# ASS_ref.fasta = reference transcriptome (Trinity assembly) used for ORF prediction.

# 1) ORF prediction
TransDecoder.LongOrfs -t ASS_ref.fasta
TransDecoder.Predict  -t ASS_ref.fasta

# 2) remove the '*' (stop codon) characters from the .pep file; InterProScan does not accept them
sed 's/\*//g' ASS_ref.fasta.transdecoder.pep > ASS_ref_clean.pep

# 3) InterProScan (limit RAM through the JVM; threads through -cpu)
export JAVA_OPTS="-Xmx16G"
cd [diretório_interproscan]
./interproscan.sh -cpu 8 \
    -i ASS_ref_clean.pep \
    -f tsv,xml,gff3 \
    -b resultados/ASS_anotacao \
    -appl cdd,hamap,ncbifam,pirsf,pfam,prints,prositepatterns,prositeprofiles,panther,smart,sfld,funfam,gene3d,superfamily \
    -goterms -pa

# flags: -goterms = Gene Ontology terms; -pa = pathway annotations (KEGG/Reactome/MetaCyc);
#        -f tsv is the table used downstream in R (column V14 = GO terms).
