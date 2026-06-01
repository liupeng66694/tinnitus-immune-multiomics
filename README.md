# tinnitus-immune-multiomics# Causal effects of immune cell phenotypes on tinnitus: a multi-omics study

**An integrated bioinformatics pipeline combining Mendelian randomization, machine learning, single‑cell transcriptomics, and drug repositioning**

---

## 📖 Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Citation](#citation)
- [License](#license)

---

## 📌 Overview

Tinnitus is a prevalent auditory symptom with limited treatment options, and its immune pathophysiology remains poorly understood. In this study, we integrated **Mendelian randomization (MR)** , transcriptomics, machine learning, single‑cell analysis, and computational drug repositioning to:

1. Identify causal associations between 731 immune cell phenotypes and tinnitus risk (20 traits nominated)
2. Select five feature genes (**SLC9A1**, **CD33**, **GABBR1**, **SH2B3**, **CRIM1**) via five machine learning algorithms (RF, SVM, LASSO, KNN, NB)
3. Construct a logistic regression‑based diagnostic nomogram (AUC = 0.906 in validation set)
4. Validate causal protective effects of SLC9A1, CRIM1 and GABBR1 via eQTL MR
5. Perform single‑cell virtual knockdown to reveal a nine‑gene B‑cell‑associated signature and NF‑κB activation
6. Screen for reversal compounds (L1000FWD) and validate TPCA‑1 binding to IKK‑β via molecular docking and 100 ns MD simulation

The complete pipeline is illustrated in **Figure 1** of the manuscript.

---

## 📂 Repository Structure

---

## 🔧 Dependencies

| Tool / Package | Version | Purpose |
|---|---|---|
| R | ≥4.2.0 | Core analysis |
| TwoSampleMR | ≥0.5.6 | MR analysis |
| Seurat | ≥5.0 | Single‑cell processing |
| Python | ≥3.9 | L1000FWD query, docking |
| AutoDock Vina | 1.2.5 | Molecular docking |
| GROMACS | 2022.2 | MD simulation |

> Full session info and package versions are available in `session_info.log`.

---

## 🚀 Quick Start

```bash
# Clone the repository
git clone https://github.com/liupeng66694/tinnitus-immune-multiomics.git
cd tinnitus-immune-multiomics

# Set up R environment (if using renv)
R -e "renv::restore()"

# Run MR analysis
Rscript R/01_MR_analysis.R

# Run all analyses sequentially
bash run_all.sh

# Run drug screening pipeline
cd drug_screening && python run_screening.py


### 📌 补充：建议添加 .zenodo.json 文件（可选但推荐）

为了让 Zenodo 自动抓取正确的元数据，可以在仓库根目录添加一个 `.zenodo.json` 文件[reference:2]：

```json
{
  "title": "Causal effects of immune cell phenotypes on tinnitus: a multi-omics study with drug repositioning",
  "creators": [
    {"name": "Jingjing Liu", "affiliation": "Bethune International Peace Hospital"},
    {"name": "Yajing Sun", "affiliation": "Bethune International Peace Hospital"},
    {"name": "Xiaoming Li", "affiliation": "Bethune International Peace Hospital"},
    {"name": "Peng Liu", "affiliation": "Bethune International Peace Hospital", "orcid": "0000-0000-0000-0000"}
  ],
  "license": "MIT",
  "description": "An integrated multi-omics pipeline combining Mendelian randomization, machine learning, single-cell transcriptomics, and computational drug repositioning to identify causal immune cell phenotypes, diagnostic biomarkers, and therapeutic candidates for tinnitus.",
  "keywords": [
    "tinnitus",
    "Mendelian randomization",
    "machine learning",
    "single-cell transcriptomics",
    "drug repositioning",
    "molecular dynamics",
    "bioinformatics"
  ],
  "access_right": "open",
  "upload_type": "software"
}
