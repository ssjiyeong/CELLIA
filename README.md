# CELLIA: Curated tissue-specific marker integration for LLM-based cell type annotation

**CELLIA** is a tissue-aware, large language model (LLM)-based workflow designed to deliver reproducible and interpretable cell type annotations for single-cell RNA-seq data.
It integrates curated tissue-specific marker gene databases with structured prompting to provide cell type predictions. 

---

## 🚀 Overview

Conventional cell type annotation relies on manual curation or reference mapping, which is often slow and inconsistent across tissues.  
**CELLIA** introduces a tissue-aware, knowledge-driven workflow that:

1. Performs differentially expressed gene (DEG) analysis.  
2. Filters candidate marker genes using curated tissue-specific marker databases.  
3. Constructs structured prompts combining tissue context and top-k markers.  
4. Queries an LLM (e.g., GPT, Claude, Gemini) for cell type annotation.  
5. Collects model outputs—predicted label, rationale, and evidence score.  
6. Stores results in AnnData/JSON/CSV format for reproducible analysis.  
7. Visualizes results through an interactive web interface (UMAP viewer, dot plots, explanation panels).

---

## 📁 Repository Structure

```
CELLIA/
│
├── cellia.py             # Core workflow
├── cellia_web.py         # Web interface
│
├── run_cellia.py         # Script to run LLM-based annotation only
├── run_cellia_web.py     # Script to run full workflow (annotation + web)
│  
├── dataset/              # Example datasets (.h5ad)
│   └── CRC.h5ad          # Example AnnData file used for testing
│ 
├── database/             # Curated tissue-specific marker gene resources
│   └── Marker_DB.csv
│ 
├── requirements.txt      # Dependencies
└── README.md             # Documentation
```

---

## ⚙️ Installation

```bash
git clone https://github.com/ssjiyeong/CELLIA.git
cd CELLIA

# ---------------------------------------------------------
# OPTION A: If you have Conda installed (recommended)
# ---------------------------------------------------------
# Create and activate a Conda environment
conda create -n cellia_env python=3.9 -y
source ~/.bashrc     # If you need
conda activate cellia_env

# ---------------------------------------------------------
# OPTION B: If you do NOT have Conda installed
# ---------------------------------------------------------
# Create and activate a Python virtual environment instead
# (Use this only when Conda is unavailable)
python -m venv .venv
source .venv/bin/activate

# ---------------------------------------------------------
# Install dependencies
# ---------------------------------------------------------
pip install -r requirements.txt

# ---------------------------------------------------------
# Run CELLIA
# ---------------------------------------------------------
export API_KEY="YOUR_API_KEY"

# Full pipeline: LLM-based annotation + interactive web visualization
python run_cellia_web.py   

# LLM-based annotation only (no web interface)
python run_cellia.py
```

> Recommended: Python ≥ 3.9 with Scanpy and Clustered AnnData

---

## 🗣️ Basic Usage

```python
from CELLIA import cellia_run
import scanpy as sc

adata = sc.read_h5ad("dataset/CRC.h5ad")

adata = cellia_run(
    adata,
    tissue_db="crc|colon",
    tissue_type="human colorectal cancer (CRC)",
    api_key="YOUR_API_KEY,
    n_top_markers=15,
    model="gpt-4.1-2025-04-14"
)
```

---

## 💻 Output Example

```json
{
  "cluster_id": "C5",
  "cell_type": "Liver Sinusoidal Endothelial Cell (LSEC)",
  "marker_explanations": {
    "LYVE1": "Highly specific marker for LSECs involved in scavenging",
    "CLEC4G": "Endocytic receptor specific to LSECs for antigen uptake"
  },
  "evidence_score": {
    "score": 0.92,
    "reason": "Multiple high-specificity markers detected"
  }
}
```

---

## 🧩 Citation

If you use **CELLIA** in your work, please cite:

> Shin, J.Y. *et al.*  
> *Curated tissue-specific marker integration enhances LLM-based cell type annotation in single-cell analysis* (2025)

---

## 📬 Contact

**Author:** Jiyeong Shin \
**GitHub:** [ssjiyeong](https://github.com/ssjiyeong) \
**Email:** sssjiyeong@gmail.com 

---

## 📜 License

