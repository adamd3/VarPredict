![VarPredict](https://github.com/adamd3/VarPredict/actions/workflows/vpd_test.yml/badge.svg)

# VarPredict
A command line tool for predicting gene expression using genotypic data.

## Input
Files required:
- `genotypes_file`: contains genotypes (may be integers or floats).
- `gene_expression_file`: normalised gene expression values.
- `metadata_file`: metadata file containing covariates to be used in models.

## Usage
Basic usage is as follows:

```
VarPredict [model] \
    -o output_directory \
    -g genotypes_file \
    -c gene_expression_file \
    -m metadata_file 
```

Available models: `random-forest`, `elastic-net`.

## Example data
Example data can be found in [this repo](https://github.com/adamd3/VarPredict_test_data).
