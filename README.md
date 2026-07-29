<p align="center">
  <img src="phagepleats/resources/data/img/logo.png" alt="PhagePleats Logo" width="300"/>
</p>

# PhagePleats

**PhagePleats** is a command-line tool for taxonomic classification of dsDNA bacteriophages (*Caudoviricetes*) using structural protein profiles.

The framework underlying PhagePleats was developed from the large-scale structural analysis described in *A novel approach to Caudoviricetes taxonomy utilising whole proteome structure-structure comparison*.

This reference taxonomy was generated from:
- 4,082 reference *Caudoviricetes* genomes (ICTV)
- 445,151 predicted protein structures (ESMFold)
- Large-scale structural clustering and phylogenetic analysis 
- RED-based clade normalization
- Yielding a revised structural taxonomy comprising 159 orders, 267 families, 502 subfamilies, and 1,189 genera across *Caudoviricetes*.

<p align="center">
  <img src="phagepleats/resources/data/img/pipeline_infographic.png" alt="pipeline-infographic" width="100%"/>
</p>

<p align="center">
  <img src="phagepleats/resources/data/img/Figure_3_RED_tree_portrait.png" alt="RED tree" width="100%"/>
</p>

Structures speak in code,<br>
Viral stories fold and flow -<br>
PhagePleats dares to know. 🦠🧬

---

## 📦 Installation (recommended: mamba)

PhagePleats is distributed as an **installable CLI tool** via a conda environment  
(required for FAISS and Foldseek compatibility).

### 1. Clone the repository

```bash
git clone https://github.com/linda5mith/PhagePleats.git
cd PhagePleats
```

### 2. Create & activate the environment
```bash
mamba env create -f environment.yml
mamba activate phagepleats
```

### 3. Running PhagePleats
```bash
phagepleats run \
  --path_to_pdbs /path/to/pdbs \
  --protein_genome_map /path/to/genome_protein_mapping.csv \
  --outdir /path/to/output_directory \
  --cores 4
```

Test pdbs and metadata found in: test_data/metadata, test_data/pdbs

| Argument         | Description                                   |
| ---------------- | --------------------------------------------- |
| `--path_to_pdbs` | Directory containing folded protein PDB files |
| `--protein_genome_map`     | Genome-to-protein mapping CSV       |
| `--outdir`       | Output directory                              |
| `--cores`        | Number of CPU cores (default: 8)              |

Test run (included example data)
```bash
phagepleats run \
  --path_to_pdbs test_data/pdbs \
  --protein_genome_map test_data/metadata/genome_protein_mapping.csv \
  --outdir /tmp/phagepleats_test \
  --cores 1
```

### 4. Outputs
🧬 taxonomy_predictions.csv

Predicted taxonomy for each input genome with confidence scores.

Includes:

Order, Family, Subfamily, Genus

Corresponding probability scores (*_prob)

🧪 taxa_novelty_summary.csv

Quantitative assessment of novelty relative to known clades, including:

Closest training genome

Shared structural protein similarity

Distance metrics

Novelty flags (e.g. Potential new genus)

### 📚 Citation

If you use PhagePleats in your work, please cite:

*A novel approach to Caudoviricetes taxonomy utilising whole proteome structure-structure comparison*