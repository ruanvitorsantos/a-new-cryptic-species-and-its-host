#!/usr/bin/env bash
# Download of Aspergillus genomes from NCBI Datasets (Section 4.2.1).
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

# Environment: micromamba env with ncbi-datasets-cli (conda-forge).

eval "$(micromamba shell hook --shell bash)"
micromamba create -n ncbi_datasets -y
micromamba activate ncbi_datasets
micromamba install -c conda-forge ncbi-datasets-cli -y

datasets download genome taxon "aspergillus niger"
unzip -o ncbi_dataset.zip -d ncbi_niger

datasets download genome taxon "aspergillus awamori"
unzip -o ncbi_dataset.zip -d ncbi_awamori

datasets download genome taxon "aspergillus welwitschiae"
unzip -o ncbi_dataset.zip -d ncbi_welwitschiae

# NOTE: re-running `datasets download` overwrites ncbi_dataset.zip, hence one output folder per taxon.
# Other taxa (remaining species/sections and outgroups): repeat with the corresponding taxon name or accession.
# The final list of 153 genomes is given in Table 1 of the manuscript.

micromamba deactivate
