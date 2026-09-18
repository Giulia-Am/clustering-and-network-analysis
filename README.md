# Computational Linear Algebra — Course Projects

This repository collects three projects developed for the **Computational Linear Algebra** course, each applying core linear algebra tools (eigendecomposition, spectral analysis, iterative methods) to real-world data problems: **dimensionality reduction**, **graph-based clustering**, and **web ranking**.

| Project | Folder | Language | Core technique |
|---|---|---|---|
| [Principal Component Analysis](#1--principal-component-analysis-pca) | `PCA/` | Python (Jupyter Notebook) | PCA, K-means clustering |
| [Spectral Clustering](#2--spectral-clustering) | `spectral_clustering/` | MATLAB | Graph Laplacians, Arnoldi eigensolver |
| [PageRank](#3--pagerank) | `pagerank/` | MATLAB | Power method, Google matrix |

---

## 1. 📊 Principal Component Analysis (PCA)

**Folder:** [`PCA/`](./PCA)

### Overview
This project applies **Principal Component Analysis (PCA)** to a real survey dataset (the *Young People Survey*), covering lifestyle, interests, and personality traits. The goal is to reduce the dimensionality of the data, interpret the resulting components, and use them to uncover meaningful groups of respondents.

### Methodology
- **Preprocessing:** encoding and cleaning of categorical and numerical survey responses.
- **Variance analysis:** computation of all principal components and their cumulative explained variance, to identify how many components are needed to retain most of the information.
- **Dimensionality reduction & interpretation:** the most informative PCs are selected, named, and interpreted based on their loadings; results are visualized through score plots.
- **Clustering:** the **K-means** algorithm is applied to the PCA-reduced data to identify groups of respondents with similar profiles.
- **Cluster evaluation:**
  - *External evaluation* — comparing clusters against meaningful demographic/personality labels.
  - *Internal evaluation* — using **silhouette scores** to assess cluster cohesion and separation.

### Contents
- `HW_PCA.ipynb` — main Jupyter notebook with full analysis, code, plots, and written commentary.

### How to run
```bash
pip install numpy pandas matplotlib scikit-learn
jupyter notebook HW_PCA.ipynb
```

---

## 2. 🕸️ Spectral Clustering

**Folder:** [`spectral_clustering/`](./spectral_clustering)

### Overview
This project implements **spectral clustering from scratch in MATLAB** and applies it to the character co-occurrence network of Victor Hugo's *Les Misérables* (D. E. Knuth's Stanford GraphBase dataset — 77 nodes, 254 weighted edges). The goal is to detect meaningful character communities using the spectral properties of the graph Laplacian.

### Methodology
- **Graph construction:** parsing of a GML file into a weighted adjacency matrix and computation of the degree matrix.
- **Graph Laplacians:** construction and comparison of three formulations:
  - Unnormalized Laplacian `L = D − A`
  - Symmetric normalized Laplacian `L_sym = I − D^(-1/2) A D^(-1/2)`
  - Random-walk Laplacian `L_rw = I − D^(-1) A`
- **Eigendecomposition:** the smallest eigenvalues/eigenvectors of the Laplacian are computed efficiently via a **shift-and-invert Arnoldi method**, validated against MATLAB's built-in `eigs`.
- **Clustering:** the spectral embedding built from the eigenvectors is normalized and passed to a custom **K-means** routine.
- **Evaluation:** cluster quality is assessed with **silhouette scores**.

### Contents
- MATLAB scripts/functions for adjacency and Laplacian construction, the spectral clustering pipeline (`spectral_clustering_m`), and the Arnoldi shift-invert eigensolver.
- Dataset file (GML format) for the *Les Misérables* character network.
- Project report (PDF) with full derivations, methodology, and results.

### How to run
Open MATLAB, set the project folder as the working directory, and run the main script (see in-folder comments/report for the exact entry-point script name and parameters).

---

## 3. 🔗 PageRank

**Folder:** [`pagerank/`](./pagerank)

### Overview
This project implements the **PageRank algorithm** from scratch in MATLAB and evaluates it on small illustrative web graphs as well as a real-world hyperlink dataset (`hollins.dat`, 6,012 pages / 23,875 links), to study how link structure and normalization affect page importance.

### Methodology
- **Adjacency & stochastic matrix:** construction of the raw link matrix and conversion into a column-stochastic matrix (`makeStochastic.m`), including proper handling of **dangling nodes**.
- **Google matrix:** construction of the damped, irreducible **Google matrix** `M = d·A + (1−d)·S` (damping factor `d = 0.85`) to guarantee a unique, well-defined ranking even for disconnected graphs.
- **Power method:** implementation of the iterative **power method** (`pagerank_power.m`) with ℓ¹-normalization and a Rayleigh-quotient convergence check.
- **Validation:** all results are cross-checked against MATLAB's `eigs` dominant-eigenvector computation (`validateWithEigs.m`).
- **Experiments:**
  1. A 4-page web, and the effect of adding a 5th page linking to Page 3 (raw adjacency matrix).
  2. The same scenario analyzed with the Google matrix.
  3. A 5-page web split into two disconnected sub-webs, showing why the Google matrix is required for well-posed ranking.
  4. Scalability test on the real-world `hollins.dat` dataset.

### Contents
- MATLAB scripts/functions: `makeStochastic.m`, `pagerank_power.m`, `validateWithEigs.m`, `dispFormattedPageRank.m`.
- `hollins.dat` — real-world hyperlink dataset.
- Project report (PDF) with full methodology, matrices, and ranked results.

### How to run
Open MATLAB, set the project folder as the working directory, and run the main script (see in-folder comments/report for the exact entry-point script name and parameters).

---


