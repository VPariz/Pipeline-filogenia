# Phylogenomics Pipeline via Compression (NCD)

[![Phylogenomics CI](https://github.com/VPariz/Pipeline-filogenia/actions/workflows/ci.yml/badge.svg)](https://github.com/VPariz/Pipeline-filogenia/actions/workflows/ci.yml)
[![Poster](https://img.shields.io/badge/Poster-X--Meeting_2026-green?style=flat-square&logo=academia&logoColor=white)](https://github.com/VPariz/Pipeline-filogenia)
[![Python](https://img.shields.io/badge/Python-3.9%2B-3776AB?style=flat-square&logo=python&logoColor=white)](https://python.org)
[![License](https://img.shields.io/badge/License-MIT-yellow?style=flat-square)](LICENSE)

**Undergraduate Research Project (PIBIC/CNPq) — Instituto Militar de Engenharia (IME) & PUC-Rio**  
*Research presented at the X-Meeting 2026 (Brazilian Association of Bioinformatics & Computational Biology - AB3C)*

---

## 🧬 Project Overview

This repository contains the source code, genomic datasets, and analytical pipelines developed for the PIBIC/CNPq research project at IME and PUC-Rio:  
**"Grammar-Based Compression for Evolutionary Distance Estimation: Mammalian Phylogenetics via NCD"** (*Pariz, Garcia & Haeusler*).

The study investigates the performance of **Normalized Compression Distance (NCD)** — grounded in Algorithmic Information Theory (Kolmogorov Complexity) — for alignment-free phylogenetic tree reconstruction across complete mammalian mitochondrial genomes (mtDNA):

$$NCD(x, y) = \frac{C(xy) - \min(C(x), C(y))}{\max(C(x), C(y))}$$

where $C(x)$ represents the compressed size of sequence $x$, and $C(xy)$ is the compressed size of the concatenated sequence $xy$.

### 🗜️ Evaluated Compression Algorithms
1. **Modern Statistical Compressors:** `zstd` (level 22), `brotli` (quality 11), and `ppmd` (order 16, 64 MB memory limit).
2. **Augmented Context-Free Grammar (ACFG) Heuristic:** Proprietary `IRR-C` algorithm featuring structural homology detection.
3. **Classical Benchmark:** `bzip2` (level -9).

### 🌳 Hierarchical Clustering Methods
- **UPGMA** (Unweighted Pair Group Method with Arithmetic Mean)
- **Complete Linkage**
- **Neighbor-Joining (NJ)**

### 📊 Topological Validation & Statistical Significance
Reconstructed phylogenetic trees are rigorously benchmarked against the **Gold-Standard Bayesian Maximum Clade Credibility (MCC) Tree (VertLife, Upham et al., 2019)** derived from 10,000 trees across 31 nuclear genes:
- **RF (Robinson-Foulds)** and **nRF (Normalized Robinson-Foulds)** distance metrics.
- **CTS (Custom Topological Score):** Clade-weighted penalty accounting for deep evolutionary divergence.
- **Monte Carlo Permutation Tests:** Empirical null model establishing statistical significance ($p < 0.005$).

---

## 👥 Authors & Advisors

- **João Victor Pariz** (Undergraduate Researcher — IME)
- **Prof. Alex de Vasconcellos Garcia, D.Sc.** (Advisor — IME)
- **Prof. Edward Hermann Haeusler, Ph.D.** (Advisor — PUC-Rio)

---

## 📂 Repository Structure

```
Pipeline-filogenia/
├── dados/genomas/            # Complete mitochondrial sequences (NCBI GenBank) in FASTA/txt
├── Resultados/
│   ├── matrizes/             # Precomputed pairwise NCD distance matrices (.txt)
│   ├── arvores/              # Reconstructed Newick trees (.nwk) and MCC reference tree
│   └── resumo_metricas.csv   # Consolidated topological benchmark metrics (RF, nRF, CTS)
├── codigos/
│   ├── pipeline_basico/
│   │   ├── compressor.py     # Pairwise genomic compression and NCD matrix computation
│   │   └── gerador_arvores.py# Hierarchical clustering and Newick tree generation
│   └── Graficos/
│       ├── main.py           # Interactive topological analysis GUI
│       └── gerar_figuras_artigo.py # Figure reproduction script for publication
└── requirements.txt
```

---

## 🚀 Reproduction & Usage

### 1. Requirements
Set up a Python 3.9+ virtual environment and install dependencies:
```bash
pip install -r requirements.txt
```

### 2. Execute Full Pipeline
To compute distance matrices and generate phylogenetic trees from raw sequences:
```bash
cd codigos/pipeline_basico/
python compressor.py
python gerador_arvores.py
```

### 3. Interactive GUI & Publication Figures
Launch the interactive tree visualization dashboard:
```bash
python codigos/Graficos/main.py
```

Export high-resolution publication-ready figures:
```bash
python codigos/Graficos/gerar_figuras_artigo.py
```

---

## 📑 Citation

If you use this pipeline, distance matrices, or methods in your research, please cite:

```bibtex
@inproceedings{pariz2026xmeeting,
  title={Grammar-Based Compression for Evolutionary Distance Estimation: Mammalian Phylogenetics via NCD},
  author={Pariz, Jo{\~a}o Victor and Garcia, Alex de Vasconcellos and Haeusler, Edward Hermann},
  booktitle={X-Meeting 2026 -- Brazilian Association of Bioinformatics and Computational Biology},
  year={2026},
  organization={AB3C}
}
```

---

## 📜 License
This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.
