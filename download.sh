#!/bin/bash

# %% This works
# dependencies: conda-forge::ncbi-datasets-cli
wget "https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/assemblies_release2_v1.0.index.csv"
index="./assemblies_release2_v1.0.index.csv"
GB_COLUMN=9
gb_accessions=( $(tail -n +2 $index | awk -F',' -v col=$GB_COLUMN '{print $col}' | head -1) )

mkdir -p raw_data
cd raw_data || exit
for accession in "${gb_accessions[@]}"; do
    echo "$accession"
    mkdir -p "$accession"
    cd "$accession" || exit
    datasets download genome accession "$accession" --filename "data.zip"
    unzip "data.zip"
    rm "data.zip"
    cd ..
done
cd ..

# ./raw_data/*/*/*/*/*_genomic.fna contains correct fasta

