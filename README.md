<p align="center">
  <img src="docs/phagepleats_logo.png" alt="PhagePleats Logo" width="300"/>
</p>

# PhagePleats 🧵🦠

**PhagePleats** is a toolkit for analyzing, clustering, and classifying phage genomes based on structural protein content. It provides a streamlined pipeline for computing protein similarity, defining novel clades, and predicting taxonomy based on shared protein profiles.

## 🚀 Features

- 🧬 Compute % shared proteins between phage genomes
- 📊 Cluster genomes into clades using UMAP + Leiden/Louvain
- 🔬 Predict taxonomic ranks (Order, Family, Subfamily, Genus)
- 🧠 Trainable model based on structural protein presence/absence
- 🧪 Benchmarking and quality control tools for model predictions

## 📦 Installation

Clone the repository and set up the environment:

```bash
git clone https://github.com/linda5mith/PhagePleats.git
cd PhagePleats
conda env create -f environment.yml
conda activate phagepleats
