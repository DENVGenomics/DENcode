# DENCode: A Probabilistic Framework for Haplotype-Informed Dengue Transmission Analysis


[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python Version](https://img.shields.io/badge/Python-≥3.8-green.svg)](https://www.python.org/)
[![Build Status](https://img.shields.io/badge/status-stable-brightgreen.svg)]()


## 1. Overview

DENCode is a novel, probabilistic model framework designed to estimate the relative strength of vector-mediated transmission linkage risk ($P_{ij}$) between pairs of Dengue virus (DENV) patients.

It resolves cryptic transmission links by integrating three critical data streams:

- Within-host viral haplotype diversity.
- Clinico-epidemiological metadata (time, location, date of onset of fever, etc.).
- Environmental covariates (ambient temperature).

The resulting $P_{ij}$ matrix forms a directed, weighted transmission network, which has been shown to be substantially more informative than networks derived solely from consensus sequences.

---

## 2. Model Formulation and Components

DENCode combines two primary kernels—one epidemiological and one genetic—into a unified, normalized probability $P_{ij}$ of transmission from source case $i$ to recipient case $j$.

### Core Equation

The transmission linkage risk $P_{ij}$ is calculated as:

$$P_{ij} = \frac{K_{ij} \cdot e^{(-\gamma G_{ij})}}{\sum_{j' \in J, j' \neq i} K_{ij'} \cdot e^{(-\gamma G_{ij'})}}$$

Where $\gamma$ is a scaling constant (configured in `parameters.py`).

### 2.1. Epidemiological Transmission Kernel ($K_{ij}$)

The Epidemiological Kernel models the sequential probabilities of vector-mediated transmission, integrating temporal plausibility, spatial proximity, and temperature-dependent mosquito life cycle factors:

$$K_{ij}=\sum_{\tau=0}^D I_{ij} (\tau) \cdot \lambda_i (\tau) \cdot S_{\text{mosq}} (\tau,T_i) \cdot P_{\text{mh}} \cdot f_{\text{spatial}} (d_{ij})$$

- $S_{\text{mosq}}$: Models mosquito survival over the Extrinsic Incubation Period (EIP), using temperature-dependent polynomial functions for mortality ($\mu(T)$) and EIP length.  
- $\lambda_i$: Incorporates human-to-mosquito infection probability based on viremia kinetics $V_{(\tau)}$.  
- $f_{\text{spatial}}$: An exponential decay function $e^{(-\delta \cdot d_{ij})}$ based on the geodesic distance $d_{ij}$ and the spatial decay parameter $\delta$.

### 2.2. Genetic Divergence Term ($G_{ij}$)

The term $e^{(-\gamma G_{ij})}$ acts as a phylogenetic penalty. $G_{ij}$ is an aggregate measure of genetic dissimilarity between the full sets of haplotypes $H_i$ and $H_j$, weighted by their measured abundance $a_k$:

$$G_{ij}=\sum_{H_k\in H_i} a_k \cdot \left(\sum_{H_l\in H_j} e^{d_{\text{adj}}^{(s)} (H_k,H_l)}\right)$$

The distance $d_{\text{adj}}^{(s)}$ is the patristic distance from a Maximum Likelihood phylogeny, which is then specifically adjusted for non-synonymous mutation hotspots to prevent spurious distance inflation.

---

## 3. Replication and Installation

### 3.1. Prerequisites

A Unix-like environment (Linux/macOS) is recommended. The following tools must be installed and accessible in the system's PATH:

- **Python**: Version 3.8 or higher  
- **IQTree2**: Required for generating the phylogenetic tree and calculating patristic distances  
- **MAFFT**: Required for sequence alignment of haplotypes  

### 3.2. Installation Steps

**Clone the Repository:**

```bash
git clone https://github.com/DENVGenomics/DENcode.git
cd DENcode
```

**Create and Activate Virtual Environment:**

```bash
python3 -m venv dencode_env
source dencode_env/bin/activate
```

**Install Python Dependencies:**

> Note: A complete list is provided in `requirements.txt`.

Example installation:
```bash
pip install numpy pandas biopython scipy
```

---

## 4. Data Structure and Input

Successful replication requires input data to be rigorously structured as follows.  
All data loading and parsing is handled by `load_data.py`.

| Directory/File | Description and Format | Used By |
|----------------|------------------------|----------|
| `Haplotypes/` | Files detailing reconstructed viral haplotypes and their within-host abundance per patient. | `compute_gij.py` |
| `Metadata/` | Clinico-Epidemiological Data (CSV/JSON): Must include Case ID, Date of Onset of Fever, Location (Postcode centroid), and environmental data (e.g., Mean Ambient Temperature). | `compute_kij.py` |
| `Hotspots/` | Text files (per serotype) listing the codon positions of non-synonymous mutation hotspots. | `compute_gij.py` |
| `Patristic_distance_calculator/` | Phylogenetic trees (Newick format) generated from input sequences via IQTree2. | `distances.py` |
| `parameters.py` | Hard-coded model constants: $\beta_h$, $P_{mh}$, spatial decay $\delta$, and functional forms for EIP and mosquito mortality $\mu(T)$. Modification of this file is essential for sensitivity analysis. | `main.py`, `model.py` |

---

## 5. Execution and Analysis

### 5.1. Primary Execution

To run the main analysis and generate the $P_{ij}$ matrix for the configured serotype(s):

```bash
python main.py
```

Outputs are written to the `Matrices/` directory.

### 5.2. Replicating Validation and Sensitivity Analysis

To perform the Monte Carlo simulations that replicate the parameter variation validation detailed in the publication:

```bash
python run_param_variation.py
```

This script automates 100+ iterations of the model, randomly sampling key epidemiological parameters within biologically informed distributions (defined in `parameters.py`), generating robust mean and credible intervals for $P_{ij}$.

### 5.3. Network Analysis

The resulting $P_{ij}$ matrix is the source for all downstream network analysis (often performed using external tools like NetworkX or R packages).  
Key steps include:

- **Thresholding:** Applying probability cut-offs (e.g., $P_{ij} \ge 0.1$ and $P_{ij} \ge 0.5$) to define network edges.  
- **Centrality Calculation:** Identifying epidemiologically significant nodes using metrics like PageRank, Betweenness Centrality, and Maximum Coreness (handled partially by `hotspots.py`).

---

## 6. License and Contact

This project is released under the **MIT License**.

