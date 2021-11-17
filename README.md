# FishAtlas: build an atlas of WGD-duplicated regions in teleost

 [![Snakemake](https://img.shields.io/badge/snakemake-≥5.5.4-brightgreen.svg)](https://snakemake.bitbucket.io) ![Snakemake-Report](https://img.shields.io/badge/snakemake-report-green.svg)

## Table of content

  - [Description](#description)
  - [Usage](#usage)
  - [Authors](#authors)
  - [License](#license)
  - [References](#references)

## Description

+Desc

+Inputs

+Outputs
    
+Figures in report

## Usage

All dependencies are listed in `envs/fish_atlas.yaml` (mainly snakemake, ete3, matpotlib and seaborn).

- Create the conda environment (this may take a few minutes):
```
conda install mamba
mamba env create -f envs/fish_atlas.yaml
```

- Activate the conda environment:
```
conda activate fish_atlas
```

- Run on toy example data:
```
snakemake --configfile config_example.yaml --cores 4
```

- Generate a snakemake report after a run:

```
snakemake --configfile config_example.yaml --report report_example.html
```

- To run on a user-defined dataset, create a new config file and format your input data following the provided example.

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
