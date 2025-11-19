# run_cellia_web_only.py

"""
Launch CELLIA web interface from pre-computed results.

This script assums that:
- AnnData already contains CELLIA annotaitons (e.g., cluster labels, cell type labels)
- LLM rationales are stored in a JSON file (e.g., cellia_output/gpt_explanations_db.json)
"""

import argparse
import scanpy as sc

from cellia_web import *

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run CELLIA web interface using pre-computed results."
    )

    parser.add_argument(
        "--adata",
        type=str,
        default="dataset/CRC.h5ad",
        help="Path to input AnnData (.h5ad) file with CELLIA results.",
    )

    parser.add_argument(
        "--rationale_json",
        type=str,
        default="cellia_output/gpt_explanations_db.json",
        help="Path to JSON file containing LLM explanations (gpt_explanations_db.json).",
    )

    parser.add_argument(
        "--port",
        type=int,
        default=8051,
        help="Port for the web app (default: 8051).",
    )

    return parser.parse_args()

def main():
    args = parse_args()

    print("Launching CELLIA web (web-only mode)...")
    print(f"  AnnData        : {args.adata}")
    print(f"  Rationale JSON : {args.rationale_json}")
    print(f"  Port           : {args.port}")


    adata = sc.read_h5ad(args.adata)

    launch_cap_style_app(
        adata=adata,
        port=args.port,
        debug=True,
        rationale_json_path=args.rationale_json,
    )


if __name__ == "__main__":
    main()