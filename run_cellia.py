# run_cellia.py

import os
import argparse
import scanpy as sc
from cellia import cellia_run

def parse_args():
    parser = argparse.ArgumentParser(description="Run CELLIA LLM-based annotation.")

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
        help="Tissue DB key (comma or pipe-separated). Default: 'crc|colon'",
    )

    parser.add_argument(
        "--tissue_type",
        type=str,
        default="human colorectal cancer (CRC)",
        help="Tissue type string for LLM prompt",
    )

    parser.add_argument(
        "--n_top_markers",
        type=int,
        default=15,
        help="Number of top markers to use per cluster",
    )

    parser.add_argument(
        "--api_key",
        type=str,
        default=os.getenv("API_KEY"),
        help="API key for LLM (default: read from environment variable API_KEY)",
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

    print("Running CELLIA...")
    print(f"Using AnnData: {args.adata}")
    print(f"Tissue DB: {args.tissue_db}")
    print(f"Tissue type: {args.tissue_type}")
    print(f"Top markers: {args.n_top_markers}")
    print(f"Model           : {args.model}")

    if args.api_key is None:
        raise RuntimeError("API_KEY is not set. Use --api_key or set environment variable API_KEY.")

    # Load input data
    adata = sc.read_h5ad("dataset/CRC.h5ad")

    # Run workflow
    result = cellia_run(
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

    print(result)

if __name__ == "__main__":
    main()