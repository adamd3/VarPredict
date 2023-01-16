![VarPredict](https://github.com/adamd3/VarPredict/actions/workflows/vpd_test.yml/badge.svg)

# VarPredict
A command line tool for predicting gene expression using genotypic data.

## Input
Files required:
- `meta_file`: metadata file (see below).
- `gpa_file`: gene presence/absence file from Panaroo (see below).
- `perc`: defines the minimum percent of strains containing a gene for i

## Usage
    VarPredict [model] \
        -o output_directory \
        -g genotypes_file \
        -c gene_expression_file \
        -c metdata_file 

Available models: `random-forest`, `elastic-net`.
Example data can be found in [this repo](https://github.com/adamd3/VarPredict_test_data).
