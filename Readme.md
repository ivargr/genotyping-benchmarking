# Benchmarking of KAGE and other genotypers

This repository contains a Snakemake-pipeline for benchmarking [KAGE](https://github.com/ivargr/kage) and other genotypers. Benchmarks can be done on both real (experimental) or simulated data. Running all the experiments will take 2-3 days using 16 CPU cores for each genotyper, as some of the genotypers require 10+ hours to run. However, running all genotypers on a small simulated dataset can be done in less than an hour.

## Reproducing the experiments in the KAGE manuscript
The branch v0.0.1 is a freeze of the code used to perform the experiments used in the KAGE manuscript. The Conda yml files in that branch will specify which versions of software were used. 

## Installation
### Step 1: Intall Snakemake and Conda
Before you start, you will need both Snakemake (to run the benchmarking pipeline) and Conda (to get all the correct dependencies. [Follow the instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) to install Snakemake if you don't have Snakemake allready.

### Step 2: Clone this repository
```bash
git clone https://github.com/ivargr/genotyping-benchmarking
```

NOTE: All dependencies will automatically be installed by Conda when you run the pipeline, except for PanGenie (which is currently not available on Conda). You will need to install PanGenie manually, and edit config.yaml to specify the installation path of PanGenie.


## Running the benchmarks
The default is to use 16 CPU cores for each method, and 40 CPU-cores to create indexes etc. If you want to change this, edit config.yaml before running.

All depenencies are handled by Conda (meaning you should use `--use-conda` with Snakemake) **except** Python dependencies, as we believe it is nice to have some control over these by installing them using your chosen Python interpreter. Thus, install Python requirements first into your chosen virtual environment:

```bash
pip install -r python_requirements.txt
```

### Run on a simulated dataset
Simply run the following. This will run all the genotypers on a small simulated dataset, specified in config.yaml and create a table with the results.

```bash
snakemake -s simulated_experiment.smk --use-conda
```

If everything goes fine, a file `figure11.html` with the following result table will be generated:

```html
+------------+---------------+------------------+-----------+-------------+----------------+---------+---------+--------------+
|            | Indels recall | Indels precision | Indels F1 | SNPs recall | SNPs precision | SNPs F1 | Runtime | Memory usage |
+------------+---------------+------------------+-----------+-------------+----------------+---------+---------+--------------+
|    KAGE    |     0.692     |      0.692       |   0.692   |    0.842    |     0.889      |  0.865  |  0 min  |     3 GB     |
|  PanGenie  |     0.769     |      0.714       |   0.741   |    0.789    |     0.714      |  0.750  |  2 min  |    48 GB     |
| Bayestyper |     0.231     |      1.000       |   0.375   |    0.158    |     1.000      |  0.273  |  1 min  |     3 GB     |
|   Malva    |     0.462     |      0.500       |   0.480   |    0.684    |     0.867      |  0.765  |  1 min  |    42 GB     |
| Graphtyper |     0.077     |      1.000       |   0.143   |    0.158    |     0.167      |  0.162  |  0 min  |     0 GB     |
|    GATK    |     0.154     |      0.125       |   0.138   |    0.263    |     0.714      |  0.385  |  0 min  |     2 GB     |
+------------+---------------+------------------+-----------+-------------+----------------+---------+---------+--------------+

```


### Run all experiments
Note: This will take several days and require a lot of RAM. It is possible to pick a subset of methods by editing figures.smk.
```bash
snakemake
```