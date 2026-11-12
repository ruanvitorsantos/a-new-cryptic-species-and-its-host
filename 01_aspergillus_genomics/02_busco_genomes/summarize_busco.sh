#!/usr/bin/env bash
# Summarise the BUSCO short_summary files of all genomes in one CSV (basis of Table S2 / Figure S1).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

# Run inside the directory that holds the BUSCO result folders.

echo "Genome,Complete,Complete_Single,Complete_Duplicated,Fragmented,Missing" > busco_summary_completude_genomas.csv

for dir in */; do
  genome_name=$(basename "$dir")
  summary_file="${dir}short_summary.specific.aspergillus_odb12.${genome_name}.txt"

  if [[ -f "$summary_file" ]]; then
    complete=$(grep -oP 'C:\K[0-9.]+%' "$summary_file" | head -n 1)
    single=$(grep -oP 'S:\K[0-9.]+%' "$summary_file" | head -n 1)
    duplicated=$(grep -oP 'D:\K[0-9.]+%' "$summary_file" | head -n 1)
    fragmented=$(grep -oP 'F:\K[0-9.]+%' "$summary_file" | head -n 1)
    missing=$(grep -oP 'M:\K[0-9.]+%' "$summary_file" | head -n 1)
    echo "${genome_name},${complete},${single},${duplicated},${fragmented},${missing}" >> busco_summary_completude_genomas.csv
  else
    echo "${genome_name},,,,," >> busco_summary_completude_genomas.csv
  fi
done

echo "Resultados extraídos e salvos em busco_summary_completude_genomas.csv"
