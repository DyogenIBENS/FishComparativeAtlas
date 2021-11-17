# FishAtlas: building an atlas of WGD-duplicated regions in teleost genomes

 [![Snakemake](https://img.shields.io/badge/snakemake-≥5.13-brightgreen.svg)](https://snakemake.bitbucket.io) ![Snakemake-Report](https://img.shields.io/badge/snakemake-report-green.svg)


FishAtlas is a snakemake pipeline to trace the evolution of sister duplicated chromosomes derived from whole genome duplication (WGD) in teleost genomes.

If you use FishAtlas, please cite [ref](TODO).

## Table of content

  - [Description](#description)
  - [Usage](#usage)
  - [Authors](#authors)
  - [License](#license)
  - [References](#references)


## Description

FishAtlas takes as input:
   1. ancestral chromosomes (pre-TGD) mapped on a subset of 4 teleost genomes (see the [examples](data/MacrosyntenyTGD/), taken from [Nakatani and McLysaght 2017](https://academic.oup.com/bioinformatics/article/33/14/i369/3953974)),
   2. genes coordinates files for all studied teleosts (see the [examples](data/example/genes/)),
   3. gene trees with the genes of all studied teleosts and outgroups (see the [example](data/example/SCORPiOs_ens89_corrected_forest.nhx)),
   4. the corresponding species tree (see the [example](data/example/sptree.nwk)).

The generated fish comparative atlas is provided in a tab-delimited file with 3 columns: the unique identifier of the post-duplication gene family, all extant teleost genes in the family and the predicted post-duplication ancestral chromosome (1a, 1b, 2a...).

## Usage

All dependencies are listed in `envs/fish_atlas.yaml` and include mainly python 3.6, snakemake, ete3, matpotlib and seaborn. You can install the dependencies directly with conda, as explained below, or manually install the packages listed in `envs/fish_atlas.yaml` before running FishAtlas.

### FishAtlas on example data

- Create the conda environment:
```
conda install mamba
mamba env create -f envs/fish_atlas.yaml
```

- Activate the conda environment:
```
conda activate fish_atlas
```

- Run on toy example data (10 teleost genomes, ~ 5 minutes):
```
snakemake --configfile config_example.yaml --cores 4
```

The output file `out_example/comparative_atlas.tsv` will be generated, along with figures with genomic annotations and statistics in `out_example/figures`.


- Generate a snakemake report after a run:

```
snakemake --configfile config_example.yaml --report report_example.html
```

The snakemake report `report_example.html` will be generated.


### FishAtlas on user-defined data

To run on a user-defined dataset, create a new configuration file and format your input data following the provided example.

## Authors

* [**Elise Parey**](mailto:elise.parey@bio.ens.psl.eu)
* **Alexandra Louis**
* **Hugues Roest Crollius**
* **Camille Berthelot**

## License

This code may be freely distributed and modified under the terms of the GNU General Public License version 3 (GPL v3) and the CeCILL licence version 2 of the CNRS:

- [LICENSE-GPLv3.txt](LICENSE-GPLv3.txt)
- [LICENSE-CeCILL.txt](LICENSE-CeCILLv2.txt)

## References

FishAtlas takes as input the pre-TGD ancestral chromosomes predictions from:

- [Nakatani and McLysaght 2017](https://academic.oup.com/bioinformatics/article/33/14/i369/3953974): Nakatani Y, McLysaght A. 2017. Genomes as documents of evolutionary history: a probabilistic macrosynteny model for the reconstruction of ancestral genomes. Bioinformatics 33:i369–i378.
