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
├── cellia.py               # Core workflow
├── cellia_web.py           # Web interface
│
├── run_cellia.py           # Script to run LLM-based annotation only
├── run_cellia_web.py       # Script to run full workflow (annotation + web)
├── cellia_tutorial.ipynb   # Tutorial to CELLIA workflow
│  
├── cellia_output/          # CELLIA annotation results (JSON / CSV)
│  
├── dataset/                # Example datasets (.h5ad)
│   └── CRC.h5ad            # Example AnnData file used for testing
│ 
├── database/               # Curated tissue-specific marker gene resources
│   └── Marker_DB.csv
│ 
├── requirements.txt        # Dependencies
└── README.md               # Documentation
```

---

## ⚙️ Installation

```bash
git clone https://github.com/ssjiyeong/CELLIA.git
cd CELLIA

# ---------------------------------------------------------
# Option A: If you have Conda installed (recommended)
# ---------------------------------------------------------
# Create and activate a Conda environment
conda create -n cellia_env python=3.9 -y
source ~/.bashrc     # If you need
conda activate cellia_env

# ---------------------------------------------------------
# Option B: If you do NOT have Conda installed
# ---------------------------------------------------------
# Create and activate a Python virtual environment instead
# (Use this only when Conda is unavailable)
python -m venv .venv
source .venv/bin/activate

# ---------------------------------------------------------
# Install dependencies
# ---------------------------------------------------------
pip install -r requirements.txt
```

> Recommended: Python ≥ 3.9 with Scanpy and Clustered AnnData

---

## 🗣️ Basic Usage

### ► Command Line Usage
After installing dependencies and creating an environment,  
you can run CELLIA in three ways:
```bash
export API_KEY="YOUR_API_KEY"
```
**A. Full wrokflow (LLM-based annotation + web interface)**
```bash
python run_cellia_web.py   
```
Then open the web interface in your brower:
```text
http://localhost:port
```

**B. LLM-based annotation only**
```bash
python run_cellia.py
```

**C. Interactive interface only** \
Runs the CELLIA web interface using pre-computed results.
This requires that the input AnnData and {LLM}_explanations_db.json already follow the CELLIA annotation output format.
```bash
python run_cellia_web_only.py
```

### ► Python script / Jupyter notebooks Usage

```python
from cellia import *
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

## 📘 Tutorial

A Jupyter notebook tutorial is provided in **cellia_tutorial.ipynb**. \
It shows the full CELLIA workflow with example data.

---

## Advanced usage (override defaults via CLI arguments)

You can override the default settings using command-line arguments.

**Run full workflow (annotation + web):**
```bash
python run_cellia_web.py \
  --adata dataset/YourAnnData.h5ad \
  --tissue_db "PBMC|blood" \
  --tissue_type "human PBMCs" \
  --n_top_markers 20 \
  --api_key "YOUR_API_KEY" \
  --model "claude-sonnet-4-5" \
  --port 8060 \
  --rationale_json cellia_output/gpt_explanations_db.json
```

**Run annotation only:**
```bash
python run_cellia.py \
  --adata dataset/YourAnnData.h5ad \
  --tissue_db "lung" \
  --tissue_type "human lung tissue" \
  --n_top_markers 10 \
  --api_key "YOUR_API_KEY" \
  --model "models/gemini-2.5-flash-lite" 
  ```

---

## 💻 Output Example

```json
{
    "1": {
        "cell_type": "Liver sinusoidal endothelial cell (LSEC)",
        "marker_explanations": {
            "FCN3": "Strongly and selectively expressed in LSECs; important for immune functions in the liver vasculature.",
            "CLEC4G": "Highly specific marker of human LSECs, mediating scavenger and cell adhesion roles.",
            "CLEC4M": "Also known as L-SIGN, specifically marks LSECs in the liver.",
            "LYVE1": "Expressed in liver sinusoidal endothelial cells involved in endocytic and scavenger functions.",
            "FCN2": "Associated with LSEC immune surveillance functions."
        },
        "evidence_score": 0.95,
        "evidence_reason": "Multiple highly specific and well-validated LSEC markers present, especially CLEC4G, CLEC4M, and FCN3, which are signature genes for liver sinusoidal endothelial cells."
    },
}
```

---

## 🧩 Citation

If you use **CELLIA** in your work, please cite:

---

## 📬 Contact

**Author:** Jiyeong Shin \
**GitHub:** [ssjiyeong](https://github.com/ssjiyeong) \
**Email:** sssjiyeong@gmail.com 

---

## 📜 License

