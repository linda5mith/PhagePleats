<p align="center">
  <img src="data/img/PhagePleats.png" alt="PhagePleats Logo" width="300"/>
</p>

# PhagePleats

**PhagePleats** is a tool for taxonomically classifying dsDNA bacteriophages of the order Caudoviricetes.

Each fold, a secret,<br>
A tapestry starts to form —<br>
PhagePleats dares to know. 🦠🧬🐍

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
```

## 📊 Interpreting PhagePleats Outputs

PhagePleats generates **three main result files**:

---

### 🧬 `taxonomy_predictions.csv`

This file contains the **predicted taxonomy** for each input phage, along with confidence scores.

**Key columns:**
- `Genome`: Identifier for the input genome
- `Order`, `Family`, `Subfamily`, `Genus`: Predicted taxonomic assignments
- `Order_prob`, `Family_prob`, etc.: Confidence score (0–1) for each predicted rank

---

### 🧩 `phage_closest_hit.csv`

This file reports the **closest training genome** to each input genome based on structural protein profiles.

**Key columns:**
- `input_genome`: Query genome name  
- `closest_training_genome`: Closest match in the reference database  
- `euclidean_distance`: Euclidean distance in presence/absence space  
- `%_shared_proteins`: Jaccard similarity to closest genome  
- `Order`, `Family`, etc.: Known taxonomy of the closest training genome  

---

### 🧪 `novel_taxa_summary.csv`

This file combines predictions, clade comparisons, and intra-clade statistics to assign **novelty flags**.

**Key columns:**
- `Genome`: Input genome ID  
- `closest_training_genome`: Most similar known genome  
- `euclidean_dist_to_closest_hit`: Distance to closest genome  
- `%_shared_with_predicted_*`: Mean Jaccard similarity to predicted clade  
- `eucl_dist_to_predicted_*`: Mean Euclidean distance to predicted clade  
- `*_z_shared_proteins` and `*_z_euclidean_distance`: Z-scores comparing your genome to typical members of the predicted clade  
- `*_novelty_flag`: Final novelty status for each rank:

  - `Likely member`  
  - `Potential new genus`, `Potential new family`, etc.  
  - `Unknown` (if intra-clade stats are unavailable)

---

These flags aim to **highlight potential novel taxa** and provide a quantitative basis for exploring viral novelty in your dataset.
