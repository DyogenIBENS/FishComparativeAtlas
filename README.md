# PARALOGY MAP: draw teleost genomes colored by post-TGD ancestral chromosomes

 [![Snakemake](https://img.shields.io/badge/snakemake-≥5.5.4-brightgreen.svg)](https://snakemake.bitbucket.io)

## Table of content

  - [Description](#description)
  - [Usage](#usage)
  - [References](#references)

## Description

This short pipeline takes as 4 inputs: (1) the predefined pre-TGD ancestral chromosomes from [ref], (2) the genes positions files, (3) the gene trees, (4) the species tree. PARALOGY MAP then assigns post-duplication chromosomes to genes of all extant TGD-duplicated species.

The major steps are the following:

- Define paralogous segments independently within each of the 4 reference species in [ref], to obtain post-duplication chromosomes annotations from pre-duplication chromosomes
- Define a consistent nomenclature across the 4 reference species, so that orthologous regions have same the post-duplication chromosomes name
- Annotate post-duplication ancestral genes in gene trees through a majority vote of reference species genes below the node
- Annotate all duplicated species using ancestral genes and draw the figures  


## Usage

All dependencies are listed in `envs/paralogy_map.yaml` (mainly snakemake, ete3, matpotlib and seaborn).
They can be easily installed via conda. One exception are some functions that I re-used from SCORPiOs, so SCORPiOs should be in the PYTHONPATH.

- Create the conda environment (this may take a few minutes)
```
conda env create -f envs/paralogy_map.yaml
```

- activate the conda environment
```
conda activate paralogy_map
```

- add SCORPiOs to the PYTHONPATH
```
export PYTHONPATH="$PYTHONPATH:/users/ldog/parey/ws1/Projects/SCORPiOs"
```

- To run on ensembl89
```
snakemake --configfile config_map_ens89.yaml --cores 4
```

- To run on genofishv3
```
snakemake --cores 14 --configfile config_map_genofish.yaml
```

- To run on a user-defined dataset, create a new config file following example in one of the two provided configfile

## References

PARALOGY MAP is based on pre-defined pre-TGD ancestral chromosomes from:

- [(Nakatani and McLysaght 2017)](https://academic.oup.com/bioinformatics/article/33/14/i369/3953974): Nakatani Y, McLysaght A. 2017. Genomes as documents of evolutionary history: a probabilistic macrosynteny model for the reconstruction of ancestral genomes. Bioinformatics 33:i369–i378.
