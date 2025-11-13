# run_cellia_web.py

import os
import argparse
import scanpy as sc

from cellia import *
from cellia_web import *

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run CELLIA full pipeline (LLM annotation + web visualization)."
    )

    parser.add_argument(
        "--adata",
        type=str,
        default="dataset/CRC.h5ad",
        help="Path to input AnnData (.h5ad) file",
    )

    parser.add_argument(
        "--tissue_db",
        type=str,
        default="crc|colon",
        help="Tissue DB key (e.g. 'crc|colon', 'lung', etc.)",
    )

    parser.add_argument(
        "--tissue_type",
        type=str,
        default="human colorectal cancer (CRC)",
        help="Tissue type string for LLM prompt.",
    )

    parser.add_argument(
        "--n_top_markers",
        type=int,
        default=15,
        help="Number of top markers to use per cluster.",
    )

    parser.add_argument(
        "--api_key",
        type=str,
        default=os.getenv("API_KEY"),
        help="API key for LLM (default: read from environment variable API_KEY).",
    )

    parser.add_argument(
        "--port",
        type=int,
        default=8051,
        help="Port for the web app (default: 8051).",
    )

    parser.add_argument(
        "--rationale_json",
        type=str,
        default="cellia_output/gpt_explanations_db.json",
        help="Path to JSON file containing LLM explanations/rationales.",
    )

    parser.add_argument(
        "--model",
        type=str,
        default="gpt-4.1-2025-04-14",
        help="Model name for LLM annotation (default: 'gpt-4.1-2025-04-14')",
    )

    return parser.parse_args()

def main():
    args = parse_args()

    print("Running CELLIA full pipeline (annotation + web)...")
    print(f"  AnnData path    : {args.adata}")
    print(f"  Tissue DB key   : {args.tissue_db}")
    print(f"  Tissue type     : {args.tissue_type}")
    print(f"  Top markers     : {args.n_top_markers}")
    print(f"  Model           : {args.model}")
    print(f"  Web port        : {args.port}")
    print(f"  Rationale JSON  : {args.rationale_json}")

    if args.api_key is None:
        raise RuntimeError(
            "API_KEY is not set. Use --api_key or set environment variable API_KEY."
        )

    # Load data
    adata = sc.read_h5ad(args.adata)

    # Run CELLIA pipeline (marker finding + filtering + LLM annotation)
    adata = cellia_run(
        adata=adata,
        tissue_db=args.tissue_db,
        tissue_type=args.tissue_type,
        api_key=args.api_key,
        n_top_markers=args.n_top_markers,
        model=args.model
    )

    adata.write_h5ad("dataset/CRC_cellia_annotated.h5ad")

    print("CELLIA annotation completed successfully.")
    print("Save to dataset/CRC_cellia_annotated.h5ad\n")
    print("Launching web interface...")

    rationale_path = "../cellia_output/gpt_explanations_db.json"
    
    # Launch interactive web app
    launch_cap_style_app(
        adata=adata,
        port=args.port,
        debug=True,
        rationale_json_path=args.rationale_json,
    )



if __name__ == "__main__":
    main()