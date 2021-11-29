# Benchmarking of KAGE and other genotypers

This repository contains a Snakemake-pipeline for benchmarking KAGE and other genotypers. Benchmarks can be done on both real (experimental) or simulated data. Running all the experiments will take 2-3 days using 16 CPU cores for each genotyper, as some of the genotypers require 10+ hours to run. However, running all genotypers on a small simulated dataset can be done in less than an hour.

## Installation
### Step 1: Clone this repository
```bash
git clone https://github.com/ivargr/genotyping-benchmarking
```

### Step 2: Install dependencies
There are many dependencies, and these can and should be installed with Conda:
```bash
cd genotyping-benchmarking
conda env create -f envs/env.yml
conda activate genotyping-benchmarking
```

NOTE: All dependencies can be installed with Conda, except for PanGenie. You will need to install PanGenie manually, and edit config.yaml to specify the installation path of PanGenie.


## Running the benchmarks
The default is to use 16 CPU cores for each method, and 40 CPU-cores to create indexes etc. If you want to change this, edit config.yaml.


### Run on a simulated dataset
Simply run the following:
```bash
snakemake -s simulated_experiment.smk
```

This will generate a file `figure11.html` with the following result table:

```html
+------------+---------------+------------------+-----------+-------------+----------------+---------+---------+--------------+
|            | Indels recall | Indels precision | Indels F1 | SNPs recall | SNPs precision | SNPs F1 | Runtime | Memory usage |
+------------+---------------+------------------+-----------+-------------+----------------+---------+---------+--------------+
|     Us     |     0.633     |      0.946       |   0.758   |    0.937    |     0.991      |  0.963  |  0 min  |     7 GB     |
|  Pangenie  |     0.610     |      0.964       |   0.747   |    0.903    |     0.993      |  0.946  |  3 min  |    49 GB     |
| Bayestyper |     0.568     |      0.993       |   0.723   |    0.870    |     0.997      |  0.929  |  5 min  |     3 GB     |
|    Gatk    |     0.969     |      0.827       |   0.892   |    0.988    |     0.995      |  0.991  |  3 min  |     5 GB     |
|   Malva    |     0.604     |      0.905       |   0.725   |    0.915    |     0.973      |  0.943  |  3 min  |    41 GB     |
| Graphtyper |     0.605     |      0.961       |   0.742   |    0.914    |     0.995      |  0.953  |  1 min  |     0 GB     |
+------------+---------------+------------------+-----------+-------------+----------------+---------+---------+--------------+



```


### Run all experiments
Note: This will take several days and require a lot of RAM. It is possible to pick a subset of methods by editing figures.smk.
```bash
snakemake
```