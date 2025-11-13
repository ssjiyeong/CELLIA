# run_cellia.py

import os
import scanpy as sc
from cellia import cellia_run

def main():
    print("Running CELLIA...")

    # Load input data
    adata = sc.read_h5ad("dataset/CRC.h5ad")

    # Run processing pipeline
    result = cellia_run(
        adata=adata,
        tissue_db="crc|colon",
        tissue_type="human colorectal cancer (CRC)",
        api_key=os.getenv("API_KEY"),
        n_top_markers=15
    )

    print("CELLIA completed successfully.")
    print(result)

if __name__ == "__main__":
    main()