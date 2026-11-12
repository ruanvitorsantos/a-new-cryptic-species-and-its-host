#!/usr/bin/env bash
# Environment setup for species-tree/species-delimitation analysis with BEAST2 + STACEY.
# Replace every [diretório_*] placeholder with the path used in your environment (see docs/PATHS.md).

# (Setup only; the run configuration was not recorded in the notes.)

micromamba create -n beast_env -c bioconda -c conda-forge beast2=2.6.3 -y
micromamba activate beast_env

# STACEY add-on
wget http://www.indriid.com/2023/STACEY.addon.v1.3.1.zip
unzip STACEY.addon.v1.3.1.zip -d $HOME/.beast/2.7/STACEY

# jModelTest (substitution model per locus)
wget https://github.com/ddarriba/jmodeltest2/files/665629/jmodeltest-2.1.10.tar.gz
