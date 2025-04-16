<p align="center">
  <img src="data/img/PhagePleats.png" alt="PhagePleats Logo" width="300"/>
</p>

# PhagePleats 🧵🦠

**PhagePleats** is a tool for taxonomically classifying dsDNA bacteriophages of the order Caudoviricetes.

## Features

- 🔬 Predict taxonomic ranks (Order, Family, Subfamily, Genus) of your phage
- 🧬 Compute % shared proteins and distance between input phage genomes and established taxonomic ranks

<p align="center">
  <img src="data/img/PhagePleats_Pipeline.png" alt="PhagePleats Pipeline" width="1000"/>
</p>

## 📦 Installation

Clone the repository and set up the environment:

```bash
git clone https://github.com/linda5mith/PhagePleats.git
cd PhagePleats
conda env create -f environment.yml
conda activate phagepleats
```

## ▶️ Running PhagePleats

Edit the `config.yml` file to include paths to your:

- `pdbs/` directory  
- genome-to-protein mapping file  
- desired output directory  

> **Note:** The default `config.yml` is already set up to run on the included test data.

Then run the pipeline:

```bash
snakemake --cores 4
